//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         PolygonalSection_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            18 Mar 2021
//

#if !defined(KRATOS_POLYGONAL_SECTION_H_INCLUDED )
#define  KRATOS_POLYGONAL_SECTION_H_INCLUDED

// #define USE_CGAL_FOR_TRIANGULATION
#define USE_DELAUNAY_FOR_TRIANGULATION

// System includes
#include <string>
#include <iostream>

// External includes
#if defined(USE_CGAL_FOR_TRIANGULATION) && defined(BREP_APPLICATION_USE_CGAL)
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#endif

// Project includes
#include "utilities/math_utils.h"
#include "custom_algebra/section/section.h"
#ifdef USE_DELAUNAY_FOR_TRIANGULATION
#include "custom_utilities/delaunay.h"
#endif

namespace Kratos
{
///@addtogroup PolygonalSectionApplication
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
/** Polygonal cut section
*/
class PolygonalSection : public Section
{
public:
    ///@name Type Definitions
    ///@{

    typedef Section BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::PointType PointType;
    typedef Vector VectorType;

    /// Pointer definition of PolygonalSection
    KRATOS_CLASS_POINTER_DEFINITION(PolygonalSection);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PolygonalSection() : BaseType()
    {}

    /// Copy constructor.
    PolygonalSection(PolygonalSection const& rOther)
    {}

    /// Destructor.
    virtual ~PolygonalSection()
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Clone this PolygonalSection
    virtual Section::Pointer Clone() const
    {
        return Section::Pointer(new PolygonalSection(*this));
    }

    /// Add a point to the list of vertices
    void AddVertex(const PointType& rPoint)
    {
        mVertices.push_back(rPoint);
    }

    /// Get the number of vertices
    std::size_t NumberOfVertices() const
    {
        return mVertices.size();
    }

    /// Clear the vertices
    void Clear()
    {
        mVertices.clear();
    }

    ///@}
    ///@name Access
    ///@{

    /// Compute the center of the section
    void ComputeCenter(PointType& rPoint) const final
    {
        noalias(rPoint) = ZeroVector(3);
        for (std::size_t i = 0; i < mVertices.size(); ++i)
        {
            noalias(rPoint) += mVertices[i];
        }
        rPoint /= mVertices.size();
    }

    /// Obtain a triangulation from the PolygonalSection
    int Triangulation(std::vector<PointType>& rPoints, std::vector<std::vector<std::size_t> >& Connectivities) const final
    {
        // compute a normal
        if (mVertices.size() < 3)
            // KRATOS_THROW_ERROR(std::logic_error, "Number of vertices is less than 3. This is not possible for the triangulation", "")
        {
            return -1;
        }

        VectorType V1(3), V2(3), N(3), T1(3), T2(3);
        noalias(V1) = mVertices[1] - mVertices[0];
        noalias(V2) = mVertices[2] - mVertices[0];
        // KRATOS_WATCH(V1)
        // KRATOS_WATCH(V2)

        noalias(T1) = V1 / norm_2(V1);
        noalias(N) = MathUtils<double>::CrossProduct(V1, V2);
        noalias(T2) = MathUtils<double>::CrossProduct(N, T1);
        T2 /= norm_2(T2);
        // KRATOS_WATCH(T1)
        // KRATOS_WATCH(T2)
        // KRATOS_WATCH(N)

        std::vector<double> XY;
        PointType Projection;
        double xmin = 1.0e99, xmax = -1.0e99, ymin = 1.0e99, ymax = -1.0e99;
        for (std::size_t i = 0; i < mVertices.size(); ++i)
        {
            noalias(Projection) = mVertices[i] - inner_prod(mVertices[i] - mVertices[0], N) * N;
            double x = inner_prod(Projection - mVertices[0], T1);
            double y = inner_prod(Projection - mVertices[0], T2);
            XY.push_back(x);
            XY.push_back(y);
            if (x < xmin) { xmin = x; }
            if (x > xmax) { xmax = x; }
            if (y < ymin) { ymin = y; }
            if (y > ymax) { ymax = y; }
        }

        // KRATOS_WATCH(XY.size())

#if defined(BREP_APPLICATION_USE_CGAL) && defined(USE_CGAL_FOR_TRIANGULATION)
        typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
        typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>   Vb;
        typedef CGAL::Triangulation_data_structure_2<Vb>                            Tds;
        typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                         Delaunay;
        typedef Kernel::Point_2                                                     Point2;

        std::vector< std::pair<Point2, unsigned int> > projected_points;
        for (std::size_t i = 0; i < XY.size() / 2; ++i)
        {
            projected_points.push_back( std::make_pair( Point2(XY[2 * i], XY[2 * i + 1]), i ) );
        }

        Delaunay triangulation;
        triangulation.insert(projected_points.begin(), projected_points.end());

        for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); ++fit)
        {
            Delaunay::Face_handle face = fit;
            std::vector<std::size_t> con(3);
            con[0] = face->vertex(0)->info();
            con[1] = face->vertex(1)->info();
            con[2] = face->vertex(2)->info();
            Connectivities.push_back(con);
        }

#elif defined(USE_DELAUNAY_FOR_TRIANGULATION)

        double dx = fabs(xmax - xmin);
        double dy = fabs(ymax - ymin);

        Delaunay D(xmin - dx, xmax + dx, ymin - dy, ymax + dy);
        for (std::size_t i = 0; i < XY.size() / 2; ++i)
        {
            D.addPoint(XY[2 * i], XY[2 * i + 1]);
        }
        // D.Print();

        auto triangles = D.getTriangles();
        // std::cout << "list of triangles" << std::endl;
        // KRATOS_WATCH(triangles.size())
        // for (auto t = triangles.begin(); t != triangles.end(); ++t)
        // {
        //     std::cout << t->pi.id << " " << t->pj.id << " " << t->pk.id << std::endl;
        // }

        Connectivities.resize(triangles.size());
        std::size_t cnt = 0;
        for (auto t = triangles.begin(); t != triangles.end(); ++t)
        {
            Connectivities[cnt].push_back(t->pi.id - 5);
            Connectivities[cnt].push_back(t->pj.id - 5);
            Connectivities[cnt].push_back(t->pk.id - 5);
            ++cnt;
        }

#endif

        // KRATOS_WATCH(Connectivities.size())

        rPoints.resize(mVertices.size());
        std::copy(mVertices.begin(), mVertices.end(), rPoints.begin());
        // KRATOS_WATCH(rPoints.size())

        return 0;
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
        return "PolygonalSection";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "Vertices:" << std::endl;
        for (std::size_t i = 0; i < mVertices.size(); ++i)
        {
            rOStream << " " << i << ": " << mVertices[i] << std::endl;
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

    std::vector<PointType> mVertices;

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
    PolygonalSection& operator=(PolygonalSection const& rOther);

    ///@}

}; // Class PolygonalSection

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
///@}

///@} addtogroup block

}  // namespace Kratos.

#undef USE_CGAL_FOR_TRIANGULATION
#undef USE_DELAUNAY_FOR_TRIANGULATION

#endif // KRATOS_POLYGONAL_SECTION_H_INCLUDED  defined
