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
//  Date:            24 Jun 2020
//

#if !defined(DELAUNAY_H_INCLUDED )
#define  DELAUNAY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes

namespace Kratos
{

/**
 * Delaunay triangulation
 */
class Delaunay
{
public:

    struct Point
    {
        std::size_t id;
        double x, y;

        Point()
        {}

        Point(const double& xx, const double& yy) : x(xx), y(yy)
        {}

        Point(const Point& p) : x(p.x), y(p.y), id(p.id)
        {}

        Point add(const double& a, const Point& p) const
        {
            return Point(x + a * p.x, y + a * p.y);
        }

        double distanceTo(const Point& p) const
        {
            return sqrt(pow(x - p.x, 2) + pow(y - p.y, 2));
        }

        double scalarProduct(const Point& p) const
        {
            return x * p.x + y * p.y;
        }

        double norm2() const
        {
            return sqrt(this->scalarProduct(*this));
        }

        bool operator==(const Point& p) const
        {
            return id == p.id;
        }

        bool operator<(const Point& p) const
        {
            return id < p.id;
        }

        Point& operator=(const Point& p)
        {
            x = p.x;
            y = p.y;
            id = p.id;
            return *this;
        }
    };

    struct Edge
    {
        Point p1, p2;

        Edge()
        {}

        Edge(const Point& pp1, const Point& pp2)
        {
            if (pp1 < pp2)
            {
                p1 = pp1;
                p2 = pp2;
            }
            else
            {
                p1 = pp2;
                p2 = pp1;
            }
        }

        Edge(const Edge& e) : p1(e.p1), p2(e.p2)
        {}

        bool operator==(const Edge& e) const
        {
            return ((p1 == e.p1) && (p2 == e.p2)) || ((p1 == e.p2) && (p2 == e.p1)) ;
        }

        bool operator<(const Edge& e) const
        {
            if (p1 == e.p1)
            {
                return p2 < e.p2;
            }
            else
            {
                return p1 < e.p1;
            }
        }

        Edge& operator=(const Edge& e)
        {
            p1 = e.p1;
            p2 = e.p2;
            return *this;
        }
    };

    struct Triangle
    {
        Point ccc, pi, pj, pk;
        double r;

        Triangle(const Point& ppi, const Point& ppj, const Point& ppk)
        {
            // select pi as the least index point
            if (ppi < ppj && ppi < ppk)
            {
                pi = ppi;
                pj = ppj;
                pk = ppk;
            }
            else if (ppj < ppi && ppj < ppk)
            {
                pi = ppj;
                pj = ppi;
                pk = ppk;
            }
            else if (ppk < ppi && ppk < ppj)
            {
                pi = ppk;
                pj = ppi;
                pk = ppj;
            }

            // make the triangle ccw
            double test = (pj.x - pi.x) * (pk.y - pi.y) - (pj.y - pi.y) * (pk.x - pi.x);
            if (test < 0) { std::swap(pj, pk); }

            // compute the center of circumferential circle and its radius
            Point dij = pi.add(-1.0, pj);
            Point djk = pj.add(-1.0, pk);
            double dxij = dij.x;
            double dyij = dij.y;
            double dxjk = djk.x;
            double dyjk = djk.y;
            double li2 = pi.scalarProduct(pi);
            double lj2 = pj.scalarProduct(pj);
            double lk2 = pk.scalarProduct(pk);
            double d = 2 * (dxij * dyjk - dxjk * dyij);
            double dx = (li2 - lj2) * dyjk - (lj2 - lk2) * dyij;
            double dy = (lj2 - lk2) * dxij - (li2 - lj2) * dxjk;

            ccc = Point(dx / d, dy / d);
            r = pi.distanceTo(ccc);
        }

        std::vector<Edge> getEdges() const
        {
            std::vector<Edge> edges(3);
            edges[0] = Edge(pi, pj);
            edges[1] = Edge(pj, pk);
            edges[2] = Edge(pk, pi);
            return edges;
        }

        bool hasPoint(const Point& p) const
        {
            return p == pi || p == pj || p == pk;
        }

        bool isInsideCC(const Point& p) const
        {
            return ccc.distanceTo(p) < r;
        }

        bool operator<(const Triangle& t) const
        {
            if (pi == t.pi)
            {
                if (pj == t.pj)
                {
                    return pk < t.pk;
                }
                else
                {
                    return pj < t.pj;
                }
            }
            else
            {
                return pi < t.pi;
            }
        }

        Triangle& operator=(const Triangle& t)
        {
            ccc = t.ccc;
            pi = t.pi;
            pj = t.pj;
            pk = t.pk;
            r = t.r;
            return *this;
        }
    };

    struct BoundaryCollector
    {
        std::set<Edge> edges;

        void addEdge(const Edge& e)
        {
            if (edges.find(e) != edges.end())
            {
                edges.erase(e);
            }
            else
            {
                edges.insert(e);
            }
        }

        template<class TEdgeContainerType>
        void addEdges(const TEdgeContainerType& es)
        {
            for (auto e = es.begin(); e != es.end(); ++e)
            {
                addEdge(*e);
            }
        }
    };

    /// Default constructor.
    Delaunay(const double& xmin, const double& xmax, const double& ymin, const double& ymax)
    {
        bps.resize(4);
        bps[0] = Point(xmin, ymin);
        bps[0].id = 1;
        bps[1] = Point(xmax, ymin);
        bps[1].id = 2;
        bps[2] = Point(xmax, ymax);
        bps[2].id = 3;
        bps[3] = Point(xmin, ymax);
        bps[3].id = 4;

        triangles.insert(Triangle(bps[0], bps[1], bps[2]));
        triangles.insert(Triangle(bps[0], bps[2], bps[3]));
    }

    /// Destructor.
    virtual ~Delaunay()
    {}

    void addPoint(const double& x, const double& y)
    {
        Point p(x, y);
        addPoint(p);
    }

    void addPoint(const Point& p)
    {
        if (points.find(p) == points.end())
        {
            Point pn(p);
            pn.id = points.size() + 5;

            auto at = this->getAffectedTriangles(pn);

            for (auto t = at.begin(); t != at.end(); ++t)
            {
                if (triangles.find(*t) != triangles.end())
                {
                    triangles.erase(*t);
                }
            }

            auto be = getBoundaryEdges(at);

            for (auto e = be.begin(); e != be.end(); ++e)
            {
                Triangle t(e->p1, e->p2, pn);
                triangles.insert(t);
            }

            points.insert(pn);
        }
    }

    std::vector<Triangle> getTriangles() const
    {
        std::vector<Triangle> l;

        for (auto t = triangles.begin(); t != triangles.end(); ++t)
        {
            if (!isBoundaryTriangle(*t))
            {
                l.push_back(*t);
            }
        }

        return l;
    }

    void Print() const
    {
        std::cout << "points:" << std::endl;
        for (auto p = points.begin(); p != points.end(); ++p)
        {
            std::cout << p->id - 4 << ": " << p->x << ", " << p->y << std::endl;
        }

        std::cout << "triangles:" << std::endl;
        auto tris = this->getTriangles();
        for (auto t = tris.begin(); t != tris.end(); ++t)
        {
            std::cout << t->pi.id - 4 << " " << t->pj.id - 4 << " " << t->pk.id - 4 << std::endl;
        }
    }

private:

    std::vector<Point> bps;
    std::set<Point> points;
    std::set<Triangle> triangles;

    std::vector<Triangle> getAffectedTriangles(const Point& p) const
    {
        std::vector<Triangle> a;
        for (auto t = triangles.begin(); t != triangles.end(); ++t)
        {
            if (t->isInsideCC(p))
            {
                a.push_back(*t);
            }
        }
        return a;
    }

    template<class TTriangleContainter>
    std::set<Edge> getBoundaryEdges(const TTriangleContainter& ts) const
    {
        BoundaryCollector bc;

        for (auto t = ts.begin(); t != ts.end(); ++t)
        {
            bc.addEdges(t->getEdges());
        }

        return bc.edges;
    }

    bool isBoundaryTriangle(const Triangle& t) const
    {
        for (auto p = bps.begin(); p != bps.end(); ++p)
        {
            if (t.hasPoint(*p))
            {
                return true;
            }
        }
        return false;
    }

}; // Class Delaunay

} // namespace Kratos

#endif // DELAUNAY_H_INCLUDED  defined
