//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         brep_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            25 Feb 2021
//

#if !defined(KRATOS_NODAL_LEVEL_SET_H_INCLUDED )
#define  KRATOS_NODAL_LEVEL_SET_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/progress.h"
#include "custom_algebra/level_set/interpolatory_nodal_level_set.h"
#include "custom_utilities/brep_mesh_utility.h"
#include "brep_application_variables.h"

#define ENABLE_PROFILING

#ifdef ENABLE_PROFILING
#include "utilities/openmp_utils.h"
#endif

namespace Kratos
{
///@addtogroup BRepApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** As name suggested, this level set approximates the level set value by using interpolation.
 * Hence nodal level set values shall be pre-computed at node
 * To create a nodal level set, simply create an implicit level set function and then call this level set
 * E.g:
 *  ls = CircularLevelSet(0.0, 0.0, 1.0)
 *  nls = NodalLevelSet(ls)
 */
class NodalLevelSet : public InterpolatoryNodalLevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(NodalLevelSet);

    typedef InterpolatoryNodalLevelSet BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    typedef Element::IntegrationMethod IntegrationMethod;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NodalLevelSet()
        : BaseType()
        , mpLevelSet(NULL)
        , mOpMode(0)
        , mPostfix("nodal_level_set_values")
    {}

    /// Constructor with level set
    NodalLevelSet(LevelSet::Pointer pLevelSet)
        : BaseType()
        , mpLevelSet(pLevelSet)
        , mOpMode(0)
        , mPostfix("nodal_level_set_values")
    {}

    /// Copy constructor.
    NodalLevelSet(NodalLevelSet const& rOther)
        : BaseType(rOther)
        , mpLevelSet(rOther.mpLevelSet)
        , mOpMode(rOther.mOpMode)
        , mPostfix(rOther.mPostfix)
        , mNodalNodalLevelSetValues(rOther.mNodalNodalLevelSetValues)
    {}

    /// Destructor.
    virtual ~NodalLevelSet() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Get the postfix
    std::string Postfix() const
    {
        return mPostfix;
    }

    /// Set the postfix
    void SetPostfix(const std::string& postfix)
    {
        mPostfix = postfix;
    }

    /// Get the operation mode
    int OperationMode() const
    {
        return mOpMode;
    }

    /// Set the operation mode:
    /// + 0: compute nodal level set values
    /// + 1: compute nodal level set values and store to file
    /// + 2: read nodal set values from file
    /// + 3: read nodal set values using GetSolutionStepValue
    void SetOperationMode(const int opMode)
    {
        mOpMode = opMode;
    }

    /// Clone this level set
    LevelSet::Pointer CloneLevelSet() const override
    {
        return LevelSet::Pointer(new NodalLevelSet(mpLevelSet->CloneLevelSet()));
    }

    /// Initialize the internal container to store the nodal level set values
    virtual void Initialize(const ModelPart::ElementsContainerType& rElements, const int configuration)
    {
        mNodalNodalLevelSetValues.clear();

        if ((mOpMode == 0) || (mOpMode == 1))
        {
            if (mpLevelSet == NULL)
            {
                KRATOS_ERROR << "The master level set is not set";
            }

#ifdef ENABLE_PROFILING
            double start = OpenMPUtils::GetCurrentTime();
#endif

            // std::set<std::size_t> nodes;
            // for (ModelPart::ElementsContainerType::const_iterator it = rElements.begin(); it != rElements.end(); ++it)
            // {
            //     for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
            //         nodes.insert(it->GetGeometry()[i].Id());
            // }
            // KRATOS_WATCH(nodes.size())

            // TODO can we parallelize this process?

            /** serial calculation of the nodal level set values **/
            // Kratos::progress_display progress(rElements.size());
            // for (ModelPart::ElementsContainerType::const_iterator it = rElements.begin(); it != rElements.end(); ++it)
            // {
            //     for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
            //     {
            //         if (mNodalNodalLevelSetValues.find(it->GetGeometry()[i].Id()) == mNodalNodalLevelSetValues.end())
            //         {
            //             // std::cout << "compute for node " << it->GetGeometry()[i].Id() << std::endl;
            //             mNodalNodalLevelSetValues[it->GetGeometry()[i].Id()] = mpLevelSet->GetValue(it->GetGeometry()[i].GetInitialPosition());
            //         }
            //     }
            //     ++progress;
            // }

            /** multithreaded calculation of the nodal level set values **/

            // collect all the nodes and coordinates
            std::unordered_map<std::size_t, double> X_coords;
            std::unordered_map<std::size_t, double> Y_coords;
            std::unordered_map<std::size_t, double> Z_coords;
            for (ModelPart::ElementsContainerType::const_iterator it = rElements.begin(); it != rElements.end(); ++it)
            {
                for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
                {
                    X_coords[it->GetGeometry()[i].Id()] = it->GetGeometry()[i].GetInitialPosition()[0];
                    Y_coords[it->GetGeometry()[i].Id()] = it->GetGeometry()[i].GetInitialPosition()[1];
                    Z_coords[it->GetGeometry()[i].Id()] = it->GetGeometry()[i].GetInitialPosition()[2];
                }
            }

            std::vector<std::size_t> node_ids(X_coords.size());
            std::size_t cnt = 0;
            for (auto it : X_coords)
            {
                node_ids[cnt++] = it.first;
            }

            // Kratos::progress_display progress(node_ids.size());
            // for (auto it : node_ids)
            // {
            //     mNodalNodalLevelSetValues[it] = mpLevelSet->GetValue(X_coords[it], Y_coords[it], Z_coords[it]);
            //     ++progress;
            // }

            // create a partition of the node array
            int number_of_threads = OpenMPUtils::GetNumThreads();
            std::size_t number_of_nodes = node_ids.size();
            std::cout << "Number of threads for nodal level set " << this->Name() << " calculation: " << number_of_threads << std::endl;
            std::cout << "Number of nodes for nodal level set " << this->Name() << " calculation: " << number_of_nodes << std::endl;

            boost::numeric::ublas::vector<unsigned int> node_partition;
            OpenMPUtils::CreatePartition(number_of_threads, number_of_nodes, node_partition);
            KRATOS_WATCH(node_partition)

            std::vector<std::vector<double> > nodal_values(number_of_threads);
            for (int k = 0; k < number_of_threads; ++k)
            {
                nodal_values[k].resize(node_partition[k + 1] - node_partition[k]);
            }

            // parallel calculation of nodal level set values
            Kratos::progress_display progress(node_ids.size());
            #pragma omp parallel for
            for (int k = 0; k < number_of_threads; ++k)
            {
                std::size_t node_id, cnt1 = 0;
                for (std::size_t i = node_partition[k]; i < node_partition[k + 1]; ++i)
                {
                    node_id = node_ids[i];
                    nodal_values[k][cnt1++] = mpLevelSet->GetValue(X_coords[node_id], Y_coords[node_id], Z_coords[node_id]);
                    ++progress;
                }
            }

            // populate the nodal level set value container
            for (int k = 0; k < number_of_threads; ++k)
            {
                std::size_t node_id, cnt1 = 0;
                for (std::size_t i = node_partition[k]; i < node_partition[k + 1]; ++i)
                {
                    node_id = node_ids[i];
                    mNodalNodalLevelSetValues[node_id] = nodal_values[k][cnt1++];
                }
            }

            if (mOpMode == 1)
            {
                // write the nodal level set to file
                std::ofstream DataFile;
                std::string filename = this->Name() + "_" + mPostfix + ".txt";
                std::cout << "Writing nodal level set values to " << filename << " ..." << std::endl;
                DataFile.open(filename.c_str());
                DataFile.precision(15);
                DataFile.setf(std::ios_base::scientific, std::ios_base::floatfield);
                for (auto it = mNodalNodalLevelSetValues.begin(); it != mNodalNodalLevelSetValues.end(); ++it)
                {
                    DataFile << it->first << " " << it->second << std::endl;
                }
                DataFile.close();
                std::cout << "Write nodal level set values to " << filename << " completed" << std::endl;
            }
        }
        else if (mOpMode == 2)
        {
            // read the nodal level set from file
            std::string filename = this->Name() + "_" + mPostfix + ".txt";
            std::cout << "Reading nodal level set values from " << filename << " ..." << std::endl;
            std::ifstream DataFile;
            DataFile.open(filename.c_str());
            if (!DataFile.is_open())
            {
                KRATOS_ERROR << filename << " does not exist";
            }
            std::string prev_word, word;
            int i = 0;
            while (DataFile >> word)
            {
                if (i == 0)
                {
                    prev_word = word;
                    i = 1;
                }
                else if (i == 1)
                {
                    i = 0;
                    std::size_t node_id = std::atoi(prev_word.c_str());
                    double value = std::atof(word.c_str());
                    mNodalNodalLevelSetValues[node_id] = value;
                }
            }
            DataFile.close();
            std::cout << "Read nodal level set values from " << filename << " completed, "
                      << mNodalNodalLevelSetValues.size() << " values are read" << std::endl;
        }
        else
        {
            KRATOS_ERROR << "Invalid ot not supported operation mode " << mOpMode;
        }
    }

    /// Initialize the internal container to store the nodal level set values
    void Initialize(const ModelPart::NodesContainerType& rNodes, const int configuration)
    {
#ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
#endif

        mNodalNodalLevelSetValues.clear();

        if (mOpMode == 0)
        {
            if (mpLevelSet == NULL)
            {
                KRATOS_ERROR << "The master level set is not set";
            }

            Kratos::progress_display progress(rNodes.size());
            if (configuration == 0)
            {
                for (ModelPart::NodesContainerType::const_iterator it = rNodes.begin(); it != rNodes.end(); ++it)
                {
                    mNodalNodalLevelSetValues[it->Id()] = mpLevelSet->GetValue(it->GetInitialPosition());
                    ++progress;
                }
            }
            else if (configuration == 1)
            {
                for (ModelPart::NodesContainerType::const_iterator it = rNodes.begin(); it != rNodes.end(); ++it)
                {
                    mNodalNodalLevelSetValues[it->Id()] = mpLevelSet->GetValue(it->X(), it->Y(), it->Z());
                    ++progress;
                }
            }

#ifdef ENABLE_PROFILING
            std::cout << "Initialize nodal level set from " << rNodes.size()
                      << " nodes completed, time = " << (OpenMPUtils::GetCurrentTime() - start) << "s" << std::endl;
#else
            std::cout << "Initialize nodal level set from " << rNodes.size()
                      << " nodes completed" << std::endl;
#endif
        }
        else if (mOpMode == 3)
        {
            for (ModelPart::NodesContainerType::const_iterator it = rNodes.begin(); it != rNodes.end(); ++it)
            {
                mNodalNodalLevelSetValues[it->Id()] = it->GetSolutionStepValue(LEVEL_SET_VALUE);
            }
        }
        else
        {
            KRATOS_ERROR << "Invalid ot not supported operation mode " << mOpMode;
        }
    }

    /// inherit from BRep
    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        if (mpLevelSet != NULL)
        {
            return mpLevelSet->ProjectOnSurface(P, Proj);
        }
        else
        {
            return -1;
        }
    }

    ///@}
    ///@name Access
    ///@{

    /// inherit from BRep
    std::size_t WorkingSpaceDimension() const final
    {
        return mpLevelSet->WorkingSpaceDimension();
    }

    /// inherit from BRep
    std::size_t LocalSpaceDimension() const final
    {
        return mpLevelSet->LocalSpaceDimension();
    }

    /// Get level set value at a node
    double GetValue(const NodeType& rNode) const override
    {
        auto it = mNodalNodalLevelSetValues.find(rNode.Id());

        if (it != mNodalNodalLevelSetValues.end())
        {
            return it->second;
        }
        else
        {
            std::cout << "WARNING: node " << rNode.Id() << " does not have level set value. Consider to initialize first the model_part" << std::endl;
            return 0.0;
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Nodal Level Set";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        if (mpLevelSet != NULL)
        {
            rOStream << Info() << " of ";
            mpLevelSet->PrintInfo(rOStream);
        }
        else
        {
            rOStream << Info();
        }
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        if (mpLevelSet != NULL)
        {
            mpLevelSet->PrintData(rOStream);
            rOStream << std::endl;
        }

        rOStream << "Nodal level set values:" << std::endl;
        for (auto it = mNodalNodalLevelSetValues.begin(); it != mNodalNodalLevelSetValues.end(); ++it)
        {
            rOStream << " " << it->first << ": " << it->second << std::endl;
        }
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::string mPostfix;

    // an internal container to store the level set values at nodes
    LevelSet::Pointer mpLevelSet;
    std::unordered_map<std::size_t, double> mNodalNodalLevelSetValues;

    int mOpMode;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NodalLevelSet& operator=(NodalLevelSet const& rOther);

    ///@}

}; // Class NodalLevelSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NODAL_LEVEL_SET_H_INCLUDED  defined
