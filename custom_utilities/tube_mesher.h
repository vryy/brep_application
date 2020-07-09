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
//  Date:            8 Jul 2020
//


#if !defined(KRATOS_TUBE_MESHER_H_INCLUDED )
#define  KRATOS_TUBE_MESHER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "custom_algebra/curve/curve.h"
#include "custom_algebra/level_set/distance_to_curve_level_set.h"


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
/**
 * class to generate the mesh for lining and grouting
 */
class TubeMesher
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TubeMesher
    KRATOS_CLASS_POINTER_DEFINITION(TubeMesher);

    typedef typename Element::GeometryType GeometryType;

    typedef typename GeometryType::PointType NodeType;

    typedef typename NodeType::PointType PointType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TubeMesher(const Curve::Pointer pCurve, const std::vector<double>& r_list,
        const std::vector<std::size_t>& nsamping_layers,
        const std::size_t& nsampling_axial, const std::size_t& nsampling_radial,
        const double& rotate_angle, const double& start_angle, const double& end_angle,
        const double& tmin, const double& tmax,
        const int& type, const std::size_t& last_node_id)
    {
        mlayers = nsamping_layers.size();
        msub_layers = nsamping_layers;
        mrings = nsampling_axial;
        msegments = nsampling_radial;

        if (type == 1)
            mnodes = 8;
        else if (type == 2)
            mnodes = 20;
        else if (type == 3)
            mnodes = 27;

        melement_connectivities.resize(nsamping_layers.size());
        mcondition_connectivities.resize(r_list.size());
        mpoints.clear();
        std::size_t last_id = last_node_id;

        if (type == 1) // H8
        {
            // create the first layer of points
            DistanceToCurveLevelSet ls(pCurve, r_list[0]);
            std::vector<std::vector<PointType> > first_layer_points;
            ls.GeneratePoints(first_layer_points, nsampling_axial+1, nsampling_radial, rotate_angle + start_angle, rotate_angle + end_angle, tmin, tmax);

            this->FillPoints(first_layer_points, last_id, mpoints, mcondition_connectivities[0], type);
            last_id = last_node_id + mpoints.size();

            // create the points for each layer
            for (std::size_t i = 0; i < nsamping_layers.size(); ++i)
            {
                mcondition_connectivities[i+1] = mcondition_connectivities[i];

                melement_connectivities[i].resize(nsamping_layers[i]);

                for (std::size_t j = 0; j < nsamping_layers[i]; ++j)
                {
                    double r = r_list[i] + (r_list[i+1]-r_list[i])*(j+1) / nsamping_layers[i];

                    DistanceToCurveLevelSet ls2(pCurve, r);

                    std::vector<std::vector<PointType> > second_layer_points;
                    ls2.GeneratePoints(second_layer_points, nsampling_axial+1, nsampling_radial, rotate_angle + start_angle, rotate_angle + end_angle, tmin, tmax);

                    std::vector<std::vector<std::vector<std::size_t> > > inner_condition_connectivities;
                    this->FillPoints(second_layer_points, last_id, mpoints, inner_condition_connectivities, type);
                    last_id = last_node_id + mpoints.size();

                    // combine points to make element
                    melement_connectivities[i][j].resize(nsampling_axial);
                    for (std::size_t k = 0; k < nsampling_axial; ++k)
                    {
                        melement_connectivities[i][j][k].resize(nsampling_radial);
                        for (std::size_t l = 0; l < nsampling_radial; ++l)
                        {
                            melement_connectivities[i][j][k][l].resize(8);

                            for (int d = 0; d < 4; ++d)
                            {
                                melement_connectivities[i][j][k][l][d] = mcondition_connectivities[i+1][k][l][d];
                                melement_connectivities[i][j][k][l][d+4] = inner_condition_connectivities[k][l][d];
                            }
                        }
                    }

                    mcondition_connectivities[i+1] = inner_condition_connectivities;
                }
            }
        }
        else if (type == 2 || type == 3) // H20 or H27
        {
            // create the first layer of points
            DistanceToCurveLevelSet ls(pCurve, r_list[0]);
            std::vector<std::vector<PointType> > first_layer_points;
            KRATOS_WATCH(nsampling_radial)
            KRATOS_WATCH(nsampling_axial)
            ls.GeneratePoints(first_layer_points, 2*nsampling_axial+1, 2*nsampling_radial, rotate_angle + start_angle, rotate_angle + end_angle, tmin, tmax);
            KRATOS_WATCH(first_layer_points.size())
            KRATOS_WATCH(first_layer_points[0].size())

            this->FillPoints(first_layer_points, last_id, mpoints, mcondition_connectivities[0], type);
            last_id = last_node_id + mpoints.size();

            // create the points for each layer
            for (std::size_t i = 0; i < nsamping_layers.size(); ++i)
            {
                mcondition_connectivities[i+1] = mcondition_connectivities[i];

                melement_connectivities[i].resize(nsamping_layers[i]);

                for (std::size_t j = 0; j < nsamping_layers[i]; ++j)
                {
                    double r1 = r_list[i] + (r_list[i+1]-r_list[i])*(j+0.5) / nsamping_layers[i];
                    double r2 = r_list[i] + (r_list[i+1]-r_list[i])*(j+1) / nsamping_layers[i];

                    DistanceToCurveLevelSet ls1(pCurve, r1);
                    DistanceToCurveLevelSet ls2(pCurve, r2);

                    std::vector<std::vector<PointType> > middle_layer_points;
                    std::vector<std::vector<PointType> > second_layer_points;
                    ls1.GeneratePoints(middle_layer_points, 2*nsampling_axial+1, 2*nsampling_radial, rotate_angle + start_angle, rotate_angle + end_angle, tmin, tmax);
                    ls2.GeneratePoints(second_layer_points, 2*nsampling_axial+1, 2*nsampling_radial, rotate_angle + start_angle, rotate_angle + end_angle, tmin, tmax);

                    std::vector<std::vector<std::vector<std::size_t> > > middle_condition_connectivities;
                    std::vector<std::vector<std::vector<std::size_t> > > inner_condition_connectivities;
                    this->FillPoints(middle_layer_points, last_id, mpoints, middle_condition_connectivities, type);
                    last_id = last_node_id + mpoints.size();
                    this->FillPoints(second_layer_points, last_id, mpoints, inner_condition_connectivities, type);
                    last_id = last_node_id + mpoints.size();

                    // combine points to make element
                    melement_connectivities[i][j].resize(nsampling_axial);
                    for (std::size_t k = 0; k < nsampling_axial; ++k)
                    {
                        melement_connectivities[i][j][k].resize(nsampling_radial);
                        for (std::size_t l = 0; l < nsampling_radial; ++l)
                        {
                            if (type == 2)
                                FillH20Connectivities(melement_connectivities[i][j][k][l], mcondition_connectivities[i+1][k][l],
                                    middle_condition_connectivities[k][l], inner_condition_connectivities[k][l]);
                            else if (type == 3)
                                FillH27Connectivities(melement_connectivities[i][j][k][l], mcondition_connectivities[i+1][k][l],
                                    middle_condition_connectivities[k][l], inner_condition_connectivities[k][l]);
                        }
                    }

                    mcondition_connectivities[i+1] = inner_condition_connectivities;
                }
            }
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Unknown type", type)

        // for H20 element the lonely points have to be removed
        if (type == 2)
        {
            std::set<std::size_t> all_nodes;

            // collect all connected nodes
            for (std::size_t i = 0; i < melement_connectivities.size(); ++i)
            {
                for (std::size_t j = 0; j < melement_connectivities[i].size(); ++j)
                {
                    for (std::size_t k = 0; k < melement_connectivities[i][j].size(); ++k)
                    {
                        for (std::size_t l = 0; l < melement_connectivities[i][j][k].size(); ++l)
                        {
                            for (std::size_t m = 0; m < melement_connectivities[i][j][k][l].size(); ++m)
                            {
                                all_nodes.insert(melement_connectivities[i][j][k][l][m]);
                            }
                        }
                    }
                }
            }

            // assign new id and put to new point container
            std::map<std::size_t, std::size_t> node_id_old_to_new;
            std::size_t new_id = 0;
            std::vector<PointType> new_points(all_nodes.size());
            for (auto it = all_nodes.begin(); it != all_nodes.end(); ++it)
            {
                new_points[new_id] = mpoints[*it-1];
                node_id_old_to_new[*it] = ++new_id;
            }
            mpoints = new_points;

            // reassign id for elements
            for (std::size_t i = 0; i < melement_connectivities.size(); ++i)
            {
                for (std::size_t j = 0; j < melement_connectivities[i].size(); ++j)
                {
                    for (std::size_t k = 0; k < melement_connectivities[i][j].size(); ++k)
                    {
                        for (std::size_t l = 0; l < melement_connectivities[i][j][k].size(); ++l)
                        {
                            for (std::size_t m = 0; m < melement_connectivities[i][j][k][l].size(); ++m)
                            {
                                melement_connectivities[i][j][k][l][m] = node_id_old_to_new[melement_connectivities[i][j][k][l][m]];
                            }
                        }
                    }
                }
            }

            // reassign id for conditions
            for (std::size_t i = 0; i < mcondition_connectivities.size(); ++i)
            {
                for (std::size_t j = 0; j < mcondition_connectivities[i].size(); ++j)
                {
                    for (std::size_t k = 0; k < mcondition_connectivities[i][j].size(); ++k)
                    {
                        for (std::size_t l = 0; l < mcondition_connectivities[i][j][k].size(); ++l)
                        {
                            mcondition_connectivities[i][j][k][l] = node_id_old_to_new[mcondition_connectivities[i][j][k][l]];
                        }
                    }
                }
            }
        }

        std::cout << "TubeMesher is initialized" << std::endl;
    }

    /// Destructor.
    virtual ~TubeMesher()
    {
        std::cout << "TubeMesher is destroyed" << std::endl;
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    std::size_t NumberOfLayers() const {return mlayers;}
    std::size_t NumberOfSubLayers(const std::size_t& layer) const {return msub_layers[layer];}
    std::size_t NumberOfRings() const {return mrings;}
    std::size_t NumberOfSegments() const {return msegments;}
    std::size_t NumberOfNodesPerSegment() const {return mnodes;}

    const std::vector<PointType>& GetPoints() const {return mpoints;}

    const std::vector<std::vector<std::vector<std::vector<std::vector<std::size_t> > > > >&
    GetElements() const {return melement_connectivities;}

    const std::vector<std::vector<std::vector<std::vector<std::size_t> > > >&
    GetConditions() const {return mcondition_connectivities;}

    /// On output: msegments x mnodes
    void GetSlices(std::vector<std::vector<std::size_t> >& conditions, const std::size_t& slice,
        const std::size_t& layer, const std::size_t& sub_layer) const
    {
        // KRATOS_WATCH(mlayers)
        // KRATOS_WATCH(msub_layers[layer])
        // KRATOS_WATCH(msegments)
        // KRATOS_WATCH(mnodes)

        if (layer >= mlayers)
            KRATOS_THROW_ERROR(std::logic_error, "The layer does not exist", layer)

        if (sub_layer >= msub_layers[layer])
            KRATOS_THROW_ERROR(std::logic_error, "The sub-layer does not exist", sub_layer)

        // get the element in ring "slice"
        conditions.resize(msegments);
        if (slice == 0)
        {
            // conditions.resize(melement_connectivities[layer][sub_layer][0].size());
            // for (std::size_t l = 0; l < melement_connectivities[layer][sub_layer][0].size(); ++l) // segment
            for (std::size_t l = 0; l < msegments; ++l) // segment
            {
                if (mnodes == 8)
                {
                    conditions[l].resize(4);
                    conditions[l][0] = melement_connectivities[layer][sub_layer][0][l][0];
                    conditions[l][1] = melement_connectivities[layer][sub_layer][0][l][1];
                    conditions[l][2] = melement_connectivities[layer][sub_layer][0][l][5];
                    conditions[l][3] = melement_connectivities[layer][sub_layer][0][l][4];
                }
                else if (mnodes == 20)
                {
                    conditions[l].resize(8);
                    conditions[l][0] = melement_connectivities[layer][sub_layer][0][l][0];
                    conditions[l][1] = melement_connectivities[layer][sub_layer][0][l][1];
                    conditions[l][2] = melement_connectivities[layer][sub_layer][0][l][5];
                    conditions[l][3] = melement_connectivities[layer][sub_layer][0][l][4];
                    conditions[l][4] = melement_connectivities[layer][sub_layer][0][l][8];
                    conditions[l][5] = melement_connectivities[layer][sub_layer][0][l][13];
                    conditions[l][6] = melement_connectivities[layer][sub_layer][0][l][16];
                    conditions[l][7] = melement_connectivities[layer][sub_layer][0][l][12];
                }
                else if (mnodes == 27)
                {
                    conditions[l].resize(9);
                    conditions[l][0] = melement_connectivities[layer][sub_layer][0][l][0];
                    conditions[l][1] = melement_connectivities[layer][sub_layer][0][l][1];
                    conditions[l][2] = melement_connectivities[layer][sub_layer][0][l][5];
                    conditions[l][3] = melement_connectivities[layer][sub_layer][0][l][4];
                    conditions[l][4] = melement_connectivities[layer][sub_layer][0][l][8];
                    conditions[l][5] = melement_connectivities[layer][sub_layer][0][l][13];
                    conditions[l][6] = melement_connectivities[layer][sub_layer][0][l][16];
                    conditions[l][7] = melement_connectivities[layer][sub_layer][0][l][12];
                    conditions[l][8] = melement_connectivities[layer][sub_layer][0][l][21];
                }
            }
        }
        else
        {
            // conditions.resize(melement_connectivities[layer][sub_layer][slice-1].size());
            // for (std::size_t l = 0; l < melement_connectivities[layer][sub_layer][slice-1].size(); ++l) // segment
            for (std::size_t l = 0; l < msegments; ++l) // segment
            {
                if (mnodes == 8)
                {
                    conditions[l].resize(4);
                    conditions[l][0] = melement_connectivities[layer][sub_layer][slice-1][l][2];
                    conditions[l][1] = melement_connectivities[layer][sub_layer][slice-1][l][3];
                    conditions[l][2] = melement_connectivities[layer][sub_layer][slice-1][l][7];
                    conditions[l][3] = melement_connectivities[layer][sub_layer][slice-1][l][6];
                }
                else if (mnodes == 20)
                {
                    conditions[l].resize(8);
                    conditions[l][0] = melement_connectivities[layer][sub_layer][slice-1][l][2];
                    conditions[l][1] = melement_connectivities[layer][sub_layer][slice-1][l][3];
                    conditions[l][2] = melement_connectivities[layer][sub_layer][slice-1][l][7];
                    conditions[l][3] = melement_connectivities[layer][sub_layer][slice-1][l][6];
                    conditions[l][4] = melement_connectivities[layer][sub_layer][slice-1][l][10];
                    conditions[l][5] = melement_connectivities[layer][sub_layer][slice-1][l][15];
                    conditions[l][6] = melement_connectivities[layer][sub_layer][slice-1][l][18];
                    conditions[l][7] = melement_connectivities[layer][sub_layer][slice-1][l][14];
                }
                else if (mnodes == 27)
                {
                    conditions[l].resize(9);
                    conditions[l][0] = melement_connectivities[layer][sub_layer][slice-1][l][2];
                    conditions[l][1] = melement_connectivities[layer][sub_layer][slice-1][l][3];
                    conditions[l][2] = melement_connectivities[layer][sub_layer][slice-1][l][7];
                    conditions[l][3] = melement_connectivities[layer][sub_layer][slice-1][l][6];
                    conditions[l][4] = melement_connectivities[layer][sub_layer][slice-1][l][10];
                    conditions[l][5] = melement_connectivities[layer][sub_layer][slice-1][l][15];
                    conditions[l][6] = melement_connectivities[layer][sub_layer][slice-1][l][18];
                    conditions[l][7] = melement_connectivities[layer][sub_layer][slice-1][l][14];
                    conditions[l][8] = melement_connectivities[layer][sub_layer][slice-1][l][23];
                }
            }
        }
    }

    /// On output: msub_layers[layer] x msegments x mnodes
    void GetSlices(std::vector<std::vector<std::vector<std::size_t> > >& conditions, const std::size_t& slice,
        const std::size_t& layer) const
    {
        if (layer >= mlayers)
            KRATOS_THROW_ERROR(std::logic_error, "The layer does not exist", layer)

        conditions.resize(msub_layers[layer]);

        for (std::size_t i = 0; i < msub_layers[layer]; ++i)
        {
            this->GetSlices(conditions[i], slice, layer, i);
        }
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Tube Mesher";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    static constexpr double Pi = 3.1415926535897932384626433832795028841971693;

    ///@}
    ///@name Member Variables
    ///@{

    std::size_t mlayers;
    std::vector<std::size_t> msub_layers;
    std::size_t mrings;
    std::size_t msegments;
    std::size_t mnodes;

    std::vector<PointType> mpoints;
    std::vector<std::vector<std::vector<std::vector<std::vector<std::size_t> > > > > melement_connectivities;
        // element connectivities: layer -> sub layer -> ring -> segment
    std::vector<std::vector<std::vector<std::vector<std::size_t> > > > mcondition_connectivities;
        // condition (surface) connectivities: layer-boundary -> ring -> segment

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// Function to fill the point list and create the list of condition connectivities on surface
    void FillPoints(const std::vector<std::vector<PointType> >& layer_points, const std::size_t& last_id,
        std::vector<PointType>& points, std::vector<std::vector<std::vector<std::size_t> > >& condition_connectivities,
        const int& type) const
    {
        std::size_t num_division_1 = layer_points.size() - 1;
        std::size_t num_division_2 = layer_points[0].size() - 1;
        // KRATOS_WATCH(num_division_1)
        // KRATOS_WATCH(num_division_2)

        for (std::size_t i = 0; i < num_division_1+1; ++i)
        {
            for (std::size_t j = 0; j < num_division_2+1; ++j)
            {
                points.push_back(layer_points[i][j]);
            }
        }

        std::size_t num_1, num_2;

        std::vector<std::size_t> node;
        if (type == 1)
        {
            node.resize(4);
            num_1 = num_division_1;
            num_2 = num_division_2 + 1;
        }
        else if (type == 2)
        {
            node.resize(8);
            num_1 = num_division_1/2;
            num_2 = (num_division_2+1)/2;
        }
        else if (type == 3)
        {
            node.resize(9);
            num_1 = num_division_1/2;
            num_2 = (num_division_2+1)/2;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Invalid type", type)

        condition_connectivities.resize(num_1);

        for (std::size_t i = 0; i < num_1; ++i)
        {
            condition_connectivities[i].resize(num_2);

            for (std::size_t j = 0; j < num_2; ++j)
            {
                if (type == 1)
                {
                    node[0] = last_id + i * (num_division_2 + 1) + j + 1;
                    node[2] = last_id + (i + 1) * (num_division_2 + 1) + j + 1;
                    if (j < num_division_2)
                    {
                        node[1] = last_id + i * (num_division_2 + 1) + j + 2;
                        node[3] = last_id + (i + 1) * (num_division_2 + 1) + j + 2;
                    }
                    else
                    {
                        node[1] = last_id + i * (num_division_2 + 1) + 1;
                        node[3] = last_id + (i + 1) * (num_division_2 + 1) + 1;
                    }

                    condition_connectivities[i][j].resize(4);
                    condition_connectivities[i][j][0] = node[0];
                    condition_connectivities[i][j][1] = node[1];
                    condition_connectivities[i][j][2] = node[3];
                    condition_connectivities[i][j][3] = node[2];
                }
                else if ((type == 2) || (type == 3))
                {
                    node[0] = last_id + 2*i * (num_division_2+1) + 2*j + 1;
                    node[3] = last_id + (2*i + 1) * (num_division_2+1) + 2*j + 1;
                    node[6] = last_id + (2*i + 2) * (num_division_2+1) + 2*j + 1;

                    node[1] = last_id + 2*i * (num_division_2+1) + 2*j + 2;
                    node[4] = last_id + (2*i + 1) * (num_division_2+1) + 2*j + 2;
                    node[7] = last_id + (2*i + 2) * (num_division_2+1) + 2*j + 2;

                    if (j < num_2-1)
                    {
                        node[2] = last_id + 2*i * (num_division_2+1) + 2*j + 3;
                        node[5] = last_id + (2*i + 1) * (num_division_2+1) + 2*j + 3;
                        node[8] = last_id + (2*i + 2) * (num_division_2+1) + 2*j + 3;
                    }
                    else
                    {
                        node[2] = last_id + 2*i * (num_division_2+1) + 1;
                        node[5] = last_id + (2*i + 1) * (num_division_2+1) + 1;
                        node[8] = last_id + (2*i + 2) * (num_division_2+1) + 1;
                    }

                    // std::cout << "node:";
                    // for (int k = 0; k < node.size(); ++k)
                    //     std::cout << " " << node[k];
                    // std::cout << std::endl;

                    if (type == 2)
                    {
                        condition_connectivities[i][j].resize(8);
                        condition_connectivities[i][j][0] = node[0];
                        condition_connectivities[i][j][1] = node[2];
                        condition_connectivities[i][j][2] = node[8];
                        condition_connectivities[i][j][3] = node[6];
                        condition_connectivities[i][j][4] = node[1];
                        condition_connectivities[i][j][5] = node[5];
                        condition_connectivities[i][j][6] = node[7];
                        condition_connectivities[i][j][7] = node[3];
                    }
                    else if (type == 3)
                    {
                        condition_connectivities[i][j].resize(9);
                        condition_connectivities[i][j][0] = node[0];
                        condition_connectivities[i][j][1] = node[2];
                        condition_connectivities[i][j][2] = node[8];
                        condition_connectivities[i][j][3] = node[6];
                        condition_connectivities[i][j][4] = node[1];
                        condition_connectivities[i][j][5] = node[5];
                        condition_connectivities[i][j][6] = node[7];
                        condition_connectivities[i][j][7] = node[3];
                        condition_connectivities[i][j][8] = node[4];
                    }
                }
            }
        }
    }

    /// fill the H20 elements providing the nodes at 3 layer
    /// The nodes at each layer shall follow the Q9 convention
    void FillH20Connectivities(std::vector<std::size_t>& element, const std::vector<std::size_t>& layer1,
        const std::vector<std::size_t>& layer2, const std::vector<std::size_t>& layer3) const
    {
        element.resize(20);

        element[0] = layer1[0];
        element[1] = layer1[1];
        element[2] = layer1[2];
        element[3] = layer1[3];
        element[8] = layer1[4];
        element[9] = layer1[5];
        element[10] = layer1[6];
        element[11] = layer1[7];

        element[12] = layer2[0];
        element[13] = layer2[1];
        element[14] = layer2[2];
        element[15] = layer2[3];

        element[4] = layer3[0];
        element[5] = layer3[1];
        element[6] = layer3[2];
        element[7] = layer3[3];
        element[16] = layer3[4];
        element[17] = layer3[5];
        element[18] = layer3[6];
        element[19] = layer3[7];
    }

    /// fill the H27 elements providing the nodes at 3 layer
    /// The nodes at each layer shall follow the Q9 convention
    void FillH27Connectivities(std::vector<std::size_t>& element, const std::vector<std::size_t>& layer1,
        const std::vector<std::size_t>& layer2, const std::vector<std::size_t>& layer3) const
    {
        element.resize(27);

        element[0] = layer1[0];
        element[1] = layer1[1];
        element[2] = layer1[2];
        element[3] = layer1[3];
        element[8] = layer1[4];
        element[9] = layer1[5];
        element[10] = layer1[6];
        element[11] = layer1[7];
        element[20] = layer1[8];

        element[12] = layer2[0];
        element[13] = layer2[1];
        element[14] = layer2[2];
        element[15] = layer2[3];
        element[21] = layer2[4];
        element[22] = layer2[5];
        element[23] = layer2[6];
        element[24] = layer2[7];
        element[26] = layer2[8];

        element[4] = layer3[0];
        element[5] = layer3[1];
        element[6] = layer3[2];
        element[7] = layer3[3];
        element[16] = layer3[4];
        element[17] = layer3[5];
        element[18] = layer3[6];
        element[19] = layer3[7];
        element[25] = layer3[8];
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
    TubeMesher& operator=(TubeMesher const& rOther);

    /// Copy constructor.
    TubeMesher(TubeMesher const& rOther);

    ///@}

}; // Class TubeMesher

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, TubeMesher& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const TubeMesher& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_TUBE_MESHER_H_INCLUDED  defined
