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

