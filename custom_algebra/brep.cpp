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
//  Date:            13 Mar 2017
//


// Project includes
#include "custom_algebra/brep.h"
#include "custom_utilities/brep_mesh_utility.h"


namespace Kratos
{

BRep::BRep() : mTOL(1.0e-10)
{}

BRep::~BRep()
{}

int BRep::CutStatus(Element::Pointer p_elem, const int& configuration) const
{
    return this->CutStatus(p_elem->GetGeometry(), configuration);
}

int BRep::CutStatus(GeometryType::Pointer p_geom, const int& configuration) const
{
    return this->CutStatus(*p_geom, configuration);
}

int BRep::CutStatus(GeometryType& r_geom, const int& configuration) const
{
    if (configuration == 0)
    {
        std::vector<PointType> points(r_geom.size());
        for (std::size_t i = 0; i < r_geom.size(); ++i)
            noalias(points[i]) = r_geom[i].GetInitialPosition();
        return CutStatusOfPoints(points);
    }
    else if (configuration == 1)
    {
        return CutStatusOfPoints(r_geom);
        // REMARK: this will use the current position of node, e.g. in dynamics
    }
}

int BRep::CutStatus(const std::vector<PointType>& r_points) const
{
    return CutStatusOfPoints(r_points);
}

int BRep::CutStatusBySampling(Element::Pointer p_elem, const std::size_t& nsampling, const int& configuration) const
{
    return this->CutStatusBySampling(p_elem->GetGeometry(), nsampling, configuration);
}

int BRep::CutStatusBySampling(GeometryType::Pointer p_geom, const std::size_t& nsampling, const int& configuration) const
{
    return this->CutStatusBySampling(*p_geom, nsampling, configuration);
}

int BRep::CutStatusBySampling(GeometryType& r_geom, const std::size_t& nsampling, const int& configuration) const
{
    std::vector<PointType> SamplingPoints;
    if (configuration == 0)
        BRepMeshUtility::GenerateSamplingPoints<0>(SamplingPoints, r_geom, nsampling);
    else if (configuration == 1)
        BRepMeshUtility::GenerateSamplingPoints<1>(SamplingPoints, r_geom, nsampling);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown configuration", configuration)
    return this->CutStatus(SamplingPoints);
}

std::string BRep::CutStatusStr(const int& stat)
{
    if (stat == _CUT)
        return "CUT";
    else if (stat == _IN)
        return "IN";
    else if (stat == _OUT)
        return "OUT";
    else
        return "Undefined";
}

}  // namespace Kratos.

