#include "Rectangle.h"
#include <iostream>
using namespace std;
using namespace corenc;
using namespace Mesh;


CRectangleBasis2y::CRectangleBasis2y()
{
    m_number = 0;
    m_normal = Point(0, 0, 0);
    m_det = 0;
    createS();
    m_1dorder = m_order = 0;
    m_hx = 1;
    m_hy = 1;
    m_points.resize(4);
    m_hx = 0;
    m_hy = 0;
}
CRectangleBasis2y::CRectangleBasis2y(const Point& p1, const Point& p2, const Point& p3, const Point& p4, const int order)
{
    m_hx = p4.x - p1.x;
    m_hy = p4.y - p1.y;
    compD(p1, p2, p3);
    compNormal(p1, p2, p3);
    createS();
    m_order = order;
    m_1dorder = 1;
    m_points.resize(4);
    m_points[0] = p1;
    m_points[1] = p2;
    m_points[2] = p3;
    m_points[3] = p4;
    m_det = (p4.x - p1.x) * (p4.y - p1.y);
    m_number = 6;
}

CRectangleBasis2y::CRectangleBasis2y(const Point* p, const int order)
{
    m_number = 6;
    m_hx = p[3].x - p[0].x;
    m_hy = p[3].y - p[0].y;
    m_points.resize(4);
    for (size_t i = 0; i < 4; ++i)
        m_points[i] = p[i];
    compD(p[0], p[1], p[2]);
    compNormal(p[0], p[1], p[2]);
    createS();
    m_order = order;
    m_1dorder = 1;
    m_det = (p[3].x - p[0].x) * (p[3].y - p[0].y);
}

const int CRectangleBasis2y::createS()
{
    m_sp = 6;
    m_s = 6;
    return 0;
}

const int CRectangleBasis2y::GetNumberOfShapeFunctions() const
{
    return m_sp;
}

const Point CRectangleBasis2y::GetNormal() const
{
    return m_normal;
}

void CRectangleBasis2y::ReverseNormal()
{
    m_normal.x = -m_normal.x;
    m_normal.y = -m_normal.y;
    m_normal.z = -m_normal.z;
}

const double CRectangleBasis2y::m_x1(const double x) const
{
    return (m_points[3].x - x) / m_hx;
}

const double CRectangleBasis2y::m_x2(const double x) const
{
    return (x - m_points[0].x) / m_hx;
}

const double CRectangleBasis2y::m_y1(const double y) const
{
    return (m_points[3].y - y) / m_hy;
}

const double CRectangleBasis2y::m_y2(const double y) const
{
    return (y - m_points[0].y) / m_hy;
}

const double CRectangleBasis2y::m_y3(const double y) const
{
    auto t = (2. * (y - m_points[0].y) / m_hy - 1.);
    return 1. - t * t;
}

const double CRectangleBasis2y::m_dy3(const double y) const
{
    auto t = (2. * (y - m_points[0].y) / m_hy - 1.);
    return -4. * t / m_hy;
}

const double CRectangleBasis2y::GetShapeFunction(const int k, const Point &p) const
{
    switch (k)
    {
    case 0:
        return m_x1(p.x) * m_y1(p.y);
    case 1:
        return m_x2(p.x) * m_y1(p.y);
    case 2:
        return m_x1(p.x) * m_y2(p.y);
    case 3:
        return m_x2(p.x) * m_y2(p.y);
    case 4:
        return m_x1(p.x) * m_y3(p.y);
    case 5:
        return m_x2(p.x) * m_y3(p.y);
    default:
        break;
    }
    return 0;
}

const Point CRectangleBasis2y::GetGradShapeFunction(const int k, const Point & p) const
{
    Point m_dpsi;
    switch (k)
    {
    case 0:
        m_dpsi.x = -m_y1(p.y) / m_hx;
        m_dpsi.y = -m_x1(p.x) / m_hy;
        break;
    case 1:
        m_dpsi.x = m_y1(p.y) / m_hx;
        m_dpsi.y = -m_x2(p.x) / m_hy;
        break;
    case 2:
        m_dpsi.x = -m_y2(p.y) / m_hx;
        m_dpsi.y = m_x1(p.x) / m_hy;
        break;
    case 3:
        m_dpsi.x = m_y2(p.y) / m_hx;
        m_dpsi.y = m_x2(p.x) / m_hy;
        break;
    case 4:
        m_dpsi.x = -m_y3(p.y) / m_hx;
        m_dpsi.y = m_x1(p.x) * m_dy3(p.y);
        break;
    case 5:
        m_dpsi.x = m_y3(p.y) / m_hx;
        m_dpsi.y = m_x2(p.x) * m_dy3(p.y);
        break;
    default:
        break;
    }
    return m_dpsi;
}

CRectangleBasis2y::CRectangleBasis2y(const CRectangleBasis2y& t)
{
    m_normal = t.m_normal;
    m_det = t.m_det;
    m_order = t.m_order;
    m_1dorder = t.m_1dorder;
    m_number = t.m_number;
    m_s = t.m_s;
    m_sp = t.m_sp;
    m_points = t.m_points;
    m_hx = t.m_hx;
    m_hy = t.m_hy;
}

void CRectangleBasis2y::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
{
    double py, pz;
    m_det = (p2.y - p1.y)*(p3.z - p1.z) - (p3.y - p1.y)*(p2.z - p1.z);
    m_det *= m_det;
    py = (p3.x - p1.x)*(p2.z - p1.z) - (p2.x - p1.x)*(p3.z - p1.z);
    py *= py;
    pz = (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
    pz *= pz;
    m_det = sqrt(m_det + py + pz);
}

void CRectangleBasis2y::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
{
    Point p;
    p.x = p1.y*(p2.z - p3.z) + p2.y*(p3.z - p1.z) + p3.y*(p1.z - p2.z);
    p.y = p1.z*(p2.x - p3.x) + p2.z*(p3.x - p1.x) + p3.z*(p1.x - p2.x);
    p.z = p1.x*(p2.y - p3.y) + p2.x*(p3.y - p1.y) + p3.x*(p1.y - p2.y);
    double t = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    p.x = (p.x) / t;
    p.y = (p.y) / t;
    p.z = (p.z) / t;
    m_normal = p;
}

const double CRectangleBasis2y::GetValue(const Point& p) const
{
    return 0.;
}
const int CRectangleBasis2y::IncreaseOrder()
{
    createS();
    ++m_order;
    ++m_1dorder;
    //m_all.push_back(m_sp);
    return 0;
}
const double CRectangleBasis2y::GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
{
    //return 0.0;
    if (node < 4)
        return f(verts[node]);
    if (node == 4)
        return f(verts[node]) - 0.5 * (f(verts[0]) + f(verts[2]));
    return f(verts[node]) - 0.5 * (f(verts[1]) + f(verts[3]));
}


