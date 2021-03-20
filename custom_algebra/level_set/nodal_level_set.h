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
#include <unordered_map>
#include <fstream>


// External includes
#include <boost/progress.hpp>


// Project includes
#include "includes/model_part.h"
#include "custom_algebra/level_set/level_set.h"
#include "custom_utilities/grid_binning.h"
#include "custom_utilities/brep_mesh_utility.h"


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

/// Short class definition.
/** As name suggested, this level set approximates the level set value by using interpolation.
 * Hence nodal level set values shall be pre-computed at node
 * To create a nodal level set, simply create an implicit level set function and then call this level set
 * E.g:
 *  ls = CircularLevelSet(0.0, 0.0, 1.0)
 *  nls = NodalLevelSet(ls)
 */
class NodalLevelSet : public LevelSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalLevelSet
    KRATOS_CLASS_POINTER_DEFINITION(NodalLevelSet);

    typedef LevelSet BaseType;

    typedef BaseType::GeometryType GeometryType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PointType PointType;

    typedef BaseType::CoordinatesArrayType CoordinatesArrayType;

    typedef Element::IntegrationMethod IntegrationMethod;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NodalLevelSet(LevelSet::Pointer pLevelSet)
    : BaseType()
    , mpLevelSet(pLevelSet)
    , mpGridBinning(NULL)
    , mOpMode(0)
    , mPostfix("nodal_level_set_values")
    {}

    /// Copy constructor.
    NodalLevelSet(NodalLevelSet const& rOther)
    : BaseType(rOther)
    , mpLevelSet(rOther.mpLevelSet)
    , mpGridBinning(rOther.mpGridBinning)
    , mOpMode(rOther.mOpMode)
    , mPostfix(rOther.mPostfix)
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
    void SetOperationMode(const int& opMode)
    {
        mOpMode = opMode;
    }


    /// Set the grid binning for this level set
    void SetGridBinning(GridBinning::Pointer pGridBinning)
    {
        mpGridBinning = pGridBinning;
    }


    /// Clone this level set
    LevelSet::Pointer CloneLevelSet() const final
    {
        return LevelSet::Pointer(new NodalLevelSet(mpLevelSet->CloneLevelSet()));
    }


    /// A utility function to use in Python interface to tell this is a nodal level set
    bool IsNodalLevelSet() const
    {
        return true;
    }


    /// Initialize the internal container to store the nodal level set values
    void Initialize(const ModelPart::ElementsContainerType& rElements, const int& configuration)
    {
        mNodalNodalLevelSetValues.clear();

        if ((mOpMode == 0) || (mOpMode == 1))
        {
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

            boost::progress_display progress(rElements.size());
            for (ModelPart::ElementsContainerType::const_iterator it = rElements.begin(); it != rElements.end(); ++it)
            {
                for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
                {
                    if (mNodalNodalLevelSetValues.find(it->GetGeometry()[i].Id()) == mNodalNodalLevelSetValues.end())
                    {
                        // std::cout << "compute for node " << it->GetGeometry()[i].Id() << std::endl;
                        mNodalNodalLevelSetValues[it->GetGeometry()[i].Id()] = mpLevelSet->GetValue(it->GetGeometry()[i].GetInitialPosition());
                    }
                }
                ++progress;
            }

            #ifdef ENABLE_PROFILING
            std::cout << "Initialize nodal level set from " << rElements.size()
                      << " elements completed, time = " << (OpenMPUtils::GetCurrentTime() - start) << "s" << std::endl;
            #else
                      std::cout << "Initialize nodal level set from " << rElements.size()
                      << " elements completed" << std::endl;
            #endif

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
                    DataFile << it->first << " " << it->second << std::endl;
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
                KRATOS_THROW_ERROR(std::logic_error, filename, "does not exist")
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
    }


    /// Initialize the internal container to store the nodal level set values
    void Initialize(const ModelPart::NodesContainerType& rNodes, const int& configuration)
    {
        #ifdef ENABLE_PROFILING
        double start = OpenMPUtils::GetCurrentTime();
        #endif

        mNodalNodalLevelSetValues.clear();

        boost::progress_display progress(rNodes.size());
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


    /// inherit from BRep
    /// Check if a point is inside/outside of the BRep
    inline bool IsInside0(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        double phi = this->CalculateOnPoint(rGeometry, local_coords);
        return (phi < 0.0);
    }


    /// inherit from BRep
    /// Check if a point is inside/outside of the BRep
    bool IsInside1(const GeometryType& rGeometry, const CoordinatesArrayType& local_coords) const final
    {
        return IsInside0(rGeometry, local_coords);
    }


    /// inherit from BRep
    /// Check if a geometry is cut by the BRep by sampling the geometry
    int CutStatusBySampling(GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const final
    {
        std::vector<CoordinatesArrayType> SamplingLocalPoints;
        BRepMeshUtility::GenerateSamplingLocalPoints(SamplingLocalPoints, r_geom, nsampling);

        const std::vector<PointType> dummy;

        return this->CutStatus(r_geom, SamplingLocalPoints, dummy);
    }


    /// inherit from BRep
    /// Check if a geometry is cut by the level set
    int CutStatus(GeometryType& r_geom, const int& configuration) const final
    {
        std::vector<std::size_t> in_list, out_list, on_list;
        for(std::size_t v = 0; v < r_geom.size(); ++v)
        {
            double phi = this->GetValue(r_geom[v]);
            if(phi < -this->GetTolerance())
                in_list.push_back(v);
            else if(phi > this->GetTolerance())
                out_list.push_back(v);
            else
                on_list.push_back(v);
        }

        int stat;
        if(in_list.size() == 0 && out_list.size() == 0)
        {
            for(std::size_t v = 0; v < r_geom.size(); ++v)
                KRATOS_WATCH(r_geom[v])
            KRATOS_WATCH(in_list.size())
            KRATOS_WATCH(out_list.size())
            KRATOS_WATCH(on_list.size())
            KRATOS_WATCH(this->GetTolerance())
            KRATOS_THROW_ERROR(std::logic_error, "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.", "")
        }
        else
        {
            if(in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if(out_list.size() == 0)
            {
                stat = BRep::_IN;
                return stat;
            }

            stat = BRep::_CUT;
            return stat;
        }

        return -99; // can't come here. Just to silence the compiler.
    }


    /// inherit from BRep
    /// Check if a set of points is cut by the BRep
    /// The geometry and the local information of the points are also provided.
    /// This function allows the use of BRep defined on grid (nodes)
    /// 0: the cell is completely inside the domain bounded by BRep
    /// 1: completely outside
    /// -1: the cell is cut by BRep
    int CutStatus(const GeometryType& r_geom,
        const std::vector<CoordinatesArrayType>& r_local_points,
        const std::vector<PointType>& r_points) const final
    {
        std::vector<std::size_t> in_list, out_list, on_list;
        for(std::size_t v = 0; v < r_local_points.size(); ++v)
        {
            double phi = this->CalculateOnPoint(r_geom, r_local_points[v]);
            if(phi < -this->GetTolerance())
                in_list.push_back(v);
            else if(phi > this->GetTolerance())
                out_list.push_back(v);
            else
                on_list.push_back(v);
        }

        int stat;
        if(in_list.size() == 0 && out_list.size() == 0)
        {
            for(std::size_t v = 0; v < r_local_points.size(); ++v)
                KRATOS_WATCH(r_local_points[v])
            KRATOS_WATCH(in_list.size())
            KRATOS_WATCH(out_list.size())
            KRATOS_WATCH(on_list.size())
            KRATOS_WATCH(this->GetTolerance())
            KRATOS_THROW_ERROR(std::logic_error, "!!!FATAL ERROR!!!The geometry is degenerated. We won't handle it.", "")
        }
        else
        {
            if(in_list.size() == 0)
            {
                stat = BRep::_OUT;
                return stat;
            }

            if(out_list.size() == 0)
            {
                stat = BRep::_IN;
                return stat;
            }

            stat = BRep::_CUT;
            return stat;
        }

        return -99; // can't come here. Just to silence the compiler.
    }


    /// inherit from BRep
    /// projects a point on the surface of level_set
    int ProjectOnSurface(const PointType& P, PointType& Proj) const final
    {
        return mpLevelSet->ProjectOnSurface(P, Proj);
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


    /// inherit from LevelSet
    /// Get level set value at a point
    double GetValue(const PointType& P) const final
    {
        if (mpGridBinning == NULL)
        {
            this->PrintInfo(std::cout);
            KRATOS_THROW_ERROR(std::logic_error, "Grid Binning is not set. GetValue for point will not work", "")
        }

        auto it = mpGridBinning->find(P);
        if (it != mpGridBinning->end())
        {
            CoordinatesArrayType local_coords;
            it->GetGeometry().PointLocalCoordinates(local_coords, P);

            return this->CalculateOnPoint(it->GetGeometry(), local_coords);
        }
        else
        {
            std::cout << "WARNING: point " << P << " is not contained in the grid binning" << std::endl;
        }
    }


    /// Get level set value at a node
    double GetValue(const NodeType& rNode) const
    {
        auto it = mNodalNodalLevelSetValues.find(rNode.Id());

        if (it != mNodalNodalLevelSetValues.end())
        {
            return it->second;
        }
        else
        {
            std::cout << "WARNING: node " << rNode.Id() << " does not have level set value. Consider to initialize first the model_part" << std::endl;
        }
    }


    /// Get level set value at integration points of element
    void CalculateOnIntegrationPoints(Element::Pointer pElement, std::vector<double>& rValues) const
    {
        return this->CalculateOnIntegrationPoints(pElement->GetGeometry(), pElement->GetIntegrationMethod(), rValues);
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const final
    {
        return "Nodal Level Set";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info() << " of ";
        mpLevelSet->PrintInfo(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        mpLevelSet->PrintData(rOStream);
        rOStream << std::endl;

        rOStream << "Nodal level set values:" << std::endl;
        for (auto it = mNodalNodalLevelSetValues.begin(); it != mNodalNodalLevelSetValues.end(); ++it)
        {
            rOStream << it->first << ": " << it->second << std::endl;
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
    GridBinning::Pointer mpGridBinning;

    int mOpMode;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    /// Get level set value at integration points of geometry
    void CalculateOnIntegrationPoints(const GeometryType& rGeometry, const IntegrationMethod& ThisIntegrationMethod,
            std::vector<double>& rValues) const
    {
        // collect the nodal values
        std::vector<double> nodal_values(rGeometry.size());

        for (std::size_t i = 0; i < rGeometry.size(); ++i)
        {
            auto it = mNodalNodalLevelSetValues.find(rGeometry[i].Id());

            if (it == mNodalNodalLevelSetValues.end())
            {
                std::stringstream ss;
                ss << "Error on " << __FUNCTION__ << ": Nodal level set value at node " << rGeometry[i].Id() << " is not assigned";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            nodal_values[i] = it->second;
        }

        const GeometryType::IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints( ThisIntegrationMethod );

        const Matrix& Ncontainer = rGeometry.ShapeFunctionsValues( ThisIntegrationMethod );

        rValues.resize(integration_points.size());

        for (std::size_t j = 0; j < integration_points.size(); ++j)
        {
            rValues[j] = 0.0;

            for (std::size_t i = 0; i < rGeometry.size(); ++i)
                rValues[j] += Ncontainer(j, i) * nodal_values[i];
        }
    }


    /// Get level set value at point in geometry
    double CalculateOnPoint(const GeometryType& rGeometry, const CoordinatesArrayType& rLocalPoint) const
    {
        // collect the nodal values
        std::vector<double> nodal_values(rGeometry.size());

        for (std::size_t i = 0; i < rGeometry.size(); ++i)
        {
            auto it = mNodalNodalLevelSetValues.find(rGeometry[i].Id());

            if (it == mNodalNodalLevelSetValues.end())
            {
                std::stringstream ss;
                ss << "Error on " << __FUNCTION__ << ": Nodal level set value at node " << rGeometry[i].Id() << " is not assigned";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            nodal_values[i] = it->second;
        }

        Vector Ncontainer;
        rGeometry.ShapeFunctionsValues( Ncontainer, rLocalPoint );

        double result = 0.0;
        for (std::size_t i = 0; i < rGeometry.size(); ++i)
            result += Ncontainer(i) * nodal_values[i];

        return result;
    }

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


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, NodalLevelSet& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const NodalLevelSet& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    rOStream << std::endl;

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#ifdef ENABLE_PROFILING
#undef ENABLE_PROFILING
#endif

#endif // KRATOS_NODAL_LEVEL_SET_H_INCLUDED  defined
