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
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
}

int BRep::CutStatus(const std::vector<PointType>& r_points) const
{
    KRATOS_THROW_ERROR(std::logic_error, "Calling the base class", __FUNCTION__)
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
        BRepMeshUtility::GenerateSamplingPoints0(SamplingPoints, r_geom, nsampling);
    else if (configuration == 1)
        BRepMeshUtility::GenerateSamplingPoints(SamplingPoints, r_geom, nsampling);
    else
        KRATOS_THROW_ERROR(std::logic_error, "Unknown configuration", configuration)
    return this->CutStatus(SamplingPoints);
}

}  // namespace Kratos.

