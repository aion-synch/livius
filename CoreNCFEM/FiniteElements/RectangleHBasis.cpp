#include "Rectangle.h"
#include <random>
#include <iostream>
#include "../../CoreNCA/Matrix.h"
#include "../../CoreNCA/MatrixSkyline.h"
using namespace std;
using namespace corenc;
using namespace Mesh;
using namespace Algebra;


CRectangleHBasis::CRectangleHBasis()
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
    m_px = 0;
    m_py = 0;
}
CRectangleHBasis::CRectangleHBasis(const Point& p1, const Point& p2, const Point& p3, const Point& p4, const int order)
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
    m_number = 9;
    m_px = m_py = 1;
}

CRectangleHBasis::CRectangleHBasis(const Point& p1, const Point& p2, const Point& p3, const Point& p4, const int px, const int py)
{
    m_hx = p4.x - p1.x;
    m_hy = p4.y - p1.y;
    compD(p1, p2, p3);
    compNormal(p1, p2, p3);
    createS();
    m_order = (px + 1) * (py + 1);
    m_px = px;
    m_py = py;
    m_1dorder = 1;
    m_points.resize(4);
    m_points[0] = p1;
    m_points[1] = p2;
    m_points[2] = p3;
    m_points[3] = p4;
    m_det = (p4.x - p1.x) * (p4.y - p1.y);
    m_number = m_order;
}

CRectangleHBasis::CRectangleHBasis(const Point* p, const int order)
{
    m_number = 4;
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
    m_px = m_py = 1;
}

CRectangleHBasis::CRectangleHBasis(const Point* p, const int px, const int py)
{
    m_number = (px + 1) * (py + 1);
    m_hx = p[3].x - p[0].x;
    m_hy = p[3].y - p[0].y;
    m_points.resize(4);
    for (size_t i = 0; i < 4; ++i)
        m_points[i] = p[i];
    compD(p[0], p[1], p[2]);
    compNormal(p[0], p[1], p[2]);
    createS();
    m_order = m_number;
    m_1dorder = 1;
    m_det = (p[3].x - p[0].x) * (p[3].y - p[0].y);
    m_px = px;
    m_py = py;
}

const int CRectangleHBasis::SetOrder(const int px, const int py)
{
    m_px = px;
    m_py = py;
    m_order = (px + 1) * (py + 1);
    return 0;
}

const int CRectangleHBasis::createS()
{
    m_sp = 9;
    m_s = 9;
    return 0;
}

const int CRectangleHBasis::GetNumberOfShapeFunctions() const
{
    return m_order;
}

const Point CRectangleHBasis::GetNormal() const
{
    return m_normal;
}

void CRectangleHBasis::ReverseNormal()
{
    m_normal.x = -m_normal.x;
    m_normal.y = -m_normal.y;
    m_normal.z = -m_normal.z;
}

const double CRectangleHBasis::m_x1(const double x) const
{
    return (m_points[3].x - x) / m_hx;
}

const double CRectangleHBasis::m_x2(const double x) const
{
    return (x - m_points[0].x) / m_hx;
}

const double CRectangleHBasis::m_y1(const double y) const
{
    return (m_points[3].y - y) / m_hy;
}

const double CRectangleHBasis::m_y2(const double y) const
{
    return (y - m_points[0].y) / m_hy;
}

const double CRectangleHBasis::m_xi(const double xi, const int n) const
{
    if (n < 2)
        return 0.;
    if (n == 2)
        return 1 - xi * xi;
    return std::pow(xi, n - 2) * (1 - xi * xi);
}

const double CRectangleHBasis::m_dxi(const double xi, const int n) const
{
    if (n < 2)
        return 0.;
    if (n == 2)
        return -2 * xi;
    return (n - 2) * std::pow(xi, n - 3) * (1 - xi * xi) - 2 * std::pow(xi, n - 1);
}


const double CRectangleHBasis::GetShapeFunction(const int k, const Point &p) const
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
    default:
        break;
    }
    auto yxi = (2. * (p.y - m_points[0].y) / m_hy - 1.);
    auto xxi = (2. * (p.x - m_points[0].x) / m_hx - 1.);
    // dump
    //cout << "k: " << k << "\t";
    if (k < 3 + m_px)
    {
        //cout << "\t" << "< 3 + m_px\t" << 3 + m_px << endl;
        return m_xi(xxi, k - 2) * m_y1(p.y);
    }
    // pump
    if (k < 4 + 2 * (m_px - 1))
    {
        //cout << "\t" << "< 4 + 2 * (m_px - 1)\t" << 4 + 2 * (m_px - 1) << endl;
        return m_xi(xxi, k - 1 - m_px) * m_y2(p.y);
    }
    // left
    if (k < 3 + 2 * (m_px - 1) + m_py)
    {
        //cout << "\t" << "< 3 + 2 * (m_px - 1) + m_py\t" << 3 + 2 * (m_px - 1) + m_py << endl;
        return m_x1(p.x) * m_xi(yxi, k - 2 - 2 * (m_px - 1));
    }
    // right
    if (k < 4 + 2 * (m_px - 1) + 2 * (m_py - 1))
    {
        //cout << "\t" << "< 4 + 2 * (m_px - 1) + 2 * (m_py - 1)\t" << 4 + 2 * (m_px - 1) + 2 * (m_py - 1) << endl;
        return m_x2(p.x) * m_xi(yxi, k - 1 - 2 * (m_px - 1) - m_py);
    }
    // inside
    int px = (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) % (m_px - 1);
    return m_xi(xxi, px ? (px + 2) : 2) * m_xi(yxi, (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) / (m_px - 1) + 2);
}

const Point CRectangleHBasis::GetGradShapeFunction(const int k, const Point & p) const
{
    Point m_dpsi;
    switch (k)
    {
    case 0:
        m_dpsi.x = -m_y1(p.y) / m_hx;
        m_dpsi.y = -m_x1(p.x) / m_hy;
        return m_dpsi;
    case 1:
        m_dpsi.x = m_y1(p.y) / m_hx;
        m_dpsi.y = -m_x2(p.x) / m_hy;
        return m_dpsi;
    case 2:
        m_dpsi.x = -m_y2(p.y) / m_hx;
        m_dpsi.y = m_x1(p.x) / m_hy;
        return m_dpsi;
    case 3:
        m_dpsi.x = m_y2(p.y) / m_hx;
        m_dpsi.y = m_x2(p.x) / m_hy;
        return m_dpsi;
    default:
        break;
    }
    auto yxi = (2. * (p.y - m_points[0].y) / m_hy - 1.);
    auto xxi = (2. * (p.x - m_points[0].x) / m_hx - 1.);
    double dyxi = 2. / m_hy;
    double dxxi = 2. / m_hx;
    // dump
    //cout << "k: " << k;
    if (k < 3 + m_px)
    {
        //cout << "\t" << "< 3 + m_px\t" << 3 + m_px << endl;
        double dy = -1. / m_hy;
        double dx = m_dxi(xxi, k - 2);
        m_dpsi.x = dxxi * dx * m_y1(p.y);
        m_dpsi.y = m_xi(xxi, k - 2) * dy;
        return m_dpsi;
    }
    // pump
    if (k < 4 + 2 * (m_px - 1))
    {
        //cout << "\t" << "< 4 + 2 * (m_px - 1)\t" << 4 + 2 * (m_px - 1) << endl;
        double dy = 1. / m_hy;
        double dx = m_dxi(xxi, k - 1 - m_px);
        m_dpsi.x = dxxi * dx * m_y2(p.y);
        m_dpsi.y = m_xi(xxi, k - 1 - m_px) * dy;
        return m_dpsi;
    }
    // left
    if (k < 3 + 2 * (m_px - 1) + m_py)
    {
        //cout << "\t" << "< 3 + 2 * (m_px - 1) + m_py\t" << 3 + 2 * (m_px - 1) + m_py << endl;
        double dx = -1. / m_hx;
        double dy = dyxi * m_dxi(yxi, k - 2 - 2 * (m_px - 1));
        m_dpsi.x = m_xi(yxi, k - 2 - 2 * (m_px - 1)) * dx;
        m_dpsi.y = m_x1(p.x) * dy;
        return m_dpsi;
    }
    // right
    if (k < 4 + 2 * (m_px - 1) + 2 * (m_py - 1))
    {
        //cout << "\t" << "< 4 + 2 * (m_px - 1) + 2 * (m_py - 1)\t" << 4 + 2 * (m_px - 1) + 2 * (m_py - 1) << endl;
        double dx = 1. / m_hx;
        double dy = dyxi * m_dxi(yxi, k - 1 - 2 * (m_px - 1) - m_py);
        m_dpsi.x = m_xi(yxi, k - 1 - 2 * (m_px - 1) - m_py) * dx;
        m_dpsi.y = m_x2(p.x) * dy;
        return m_dpsi;
    }
    // inside
    //cout << "inside" << endl;
    int px = (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) % (m_px - 1);

    int p_y = (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) / (m_px - 1) + 2;

    double dx = m_dxi(xxi, px ? (px + 2) : 2) * dxxi;
    double dy = m_dxi(yxi, (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) / (m_px - 1) + 2) * dyxi;

    m_dpsi.x = m_xi(yxi, (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) / (m_px - 1) + 2) * dx;
    m_dpsi.y = m_xi(xxi, px ? (px + 2) : 2) * dy;

    //cout << (px ? (px + 2) : 2) << "\t";
    //cout << (k - 2 * (m_px - 1) - 2 * (m_py - 1) - 4) / (m_px - 1) + 2 << endl;

    return m_dpsi;
}

CRectangleHBasis::CRectangleHBasis(const CRectangleHBasis& t)
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
    m_px = t.m_px;
    m_py = t.m_py;
}

void CRectangleHBasis::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
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

void CRectangleHBasis::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
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

const double CRectangleHBasis::GetValue(const Point& p) const
{
    return 0.;
}
const int CRectangleHBasis::IncreaseOrder()
{
    createS();
    ++m_order;
    ++m_1dorder;
    //m_all.push_back(m_sp);
    return 0;
}
const double CRectangleHBasis::GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
{
    //return 0.0;
    //if (node < 4)
    //    return f(verts[node]);
    //return 0;
    /*switch (node)
    {
        case 4:
            return f(verts[node]) - 0.5 * (f(verts[0]) + f(verts[1]));
        case 5:
            return f(verts[node]) - 0.5 * (f(verts[2]) + f(verts[3]));
        case 6:
            return f(verts[node]) - 0.5 * (f(verts[0]) + f(verts[2]));
        case 7:
            return f(verts[node]) - 0.5 * (f(verts[1]) + f(verts[3]));
    }*/
    Matrix matrix;
    auto total = (m_px + 1) * (m_py + 1);
    //cout << m_px << endl;
    //cout << m_py << endl;
    //cout << total << endl;
    matrix.Create(total);
    vector<double> qq2(total);

    Point temp;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0., 1.);
    // left
    for (int i = 0; i < m_py - 1; ++i)
    {
        temp = Point(verts[0].x, verts[0].y + (verts[3].y - verts[0].y) / (i + 2));
        for (int k = 0; k < total; ++k)
        {
            matrix.AddElement(i + 4, k, GetShapeFunction(k, temp));
            //cout << GetShapeFunction(k, temp) << "\t";
        }
        qq2[i + 4] = f(temp);
    }
    // right
    for (int i = 0; i < m_py - 1; ++i)
    {
        temp = Point(verts[3].x, verts[0].y + (verts[3].y - verts[0].y) / (i + 2));
        for (int k = 0; k < total; ++k)
        {
            matrix.AddElement(m_py + i + 3, k, GetShapeFunction(k, temp));
            //cout << GetShapeFunction(k, temp) << "\t";
        }
        qq2[m_py + i + 3] = f(temp);
    }
    // dump
    for (int i = 0; i < m_px - 1; ++i)
    {
        temp = Point(verts[0].x + (verts[3].x - verts[0].x) / (i + 2), verts[0].y);
        for (int k = 0; k < total; ++k)
        {
            matrix.AddElement(2 * m_py + i + 2, k, GetShapeFunction(k, temp));
            //cout << GetShapeFunction(k, temp) << "\t";
        }
        qq2[2 * m_py + i + 2] = f(temp);
    }
    // pump
    for (int i = 0; i < m_px - 1; ++i)
    {
        temp = Point(verts[0].x + (verts[3].x - verts[0].x) / (i + 2), verts[3].y);
        for (int k = 0; k < total; ++k)
        {
            matrix.AddElement(m_px + 2 * m_py + i + 1, k, GetShapeFunction(k, temp));
            //cout << GetShapeFunction(k, temp) << "\t";
        }
        qq2[m_px + 2 * m_py + i + 1] = f(temp);
    }

    int cent = 2 * m_px + 2 * m_py;

    for (int i = cent; i < total; ++i)
    {
        //temp = Point(verts[0].x + (verts[3].x - verts[0].x) / (i + 1), verts[0].y + (verts[3].y - verts[0].y) / (i + 3));
        temp = Point(verts[0].x + dist(mt) * (verts[3].x - verts[0].x), verts[0].y + dist(mt) * (verts[3].y - verts[0].y));
        for (int j = 0; j < total; ++j)
        {
            matrix.AddElement(i, j, GetShapeFunction(j, temp));
            //cout << GetShapeFunction(j, temp) << "\t";
        }
        //cout << endl;
        qq2[i] = f(temp);
    }
    /*for (int i = 0; i < m_py + 1; ++i)
    {
        for (int j = 0; j < m_px + 1; ++j)
        {
            temp = Point(verts[0].x + (verts[3].x - verts[0].x) / (j + 1), verts[0].y + (verts[3].y - verts[0].y) / (i + 1));
            for (int k = 0; k < total; ++k)
            {
                matrix.AddElement(i * (m_px + 1) + j, k, GetShapeFunction(k, temp));
                cout << GetShapeFunction(k, temp) << "\t";
            }
            cout << endl;
            //cout << i * (m_px + 1) + j << endl;
            qq2[i * (m_px + 1) + j] = f(temp);
        }
    }*/
    /*for (int i = 0; i < total; ++i)
    {
        //temp = Point(verts[0].x + (verts[3].x - verts[0].x) / (i + 1), verts[0].y + (verts[3].y - verts[0].y) / (i + 3));
        temp = Point(verts[0].x + dist(mt) * (verts[3].x - verts[0].x), verts[0].y + dist(mt) * (verts[3].y - verts[0].y));
        for (int j = 0; j < total; ++j)
        {
            matrix.AddElement(i, j, GetShapeFunction(j, temp));
            cout << GetShapeFunction(j, temp) << "\t";
        }
        cout << endl;
        qq2[i] = f(temp);
    }*/
    //return 0.;
    qq2[0] = f(verts[0]);
    qq2[1] = f(verts[1]);
    qq2[2] = f(verts[2]);
    qq2[3] = f(verts[3]);
    //qq2[4] = 0.5 * (f(verts[0]) + f(verts[1]));
    matrix.NullRow(0);
    matrix.NullRow(1);
    matrix.NullRow(2);
    matrix.NullRow(3);
    //matrix.NullRow(4);
    /*for (int i = 0; i < total; ++i)
    {
        for (int j = 0; j < total; ++j)
        {
            cout << matrix.GetElement(i, j) << "\t";
        }
        cout << endl;
    }*/
    ESolver solver;
    vector<double> sol(total);
    solver.Gauss(matrix, qq2, sol);
    /*for (int i = 0; i < total; ++i)
        cout << qq2[i] << endl;
    //Point tp{verts[0].x + (verts[3].x - verts[0].x) / 5, verts[0].y + (verts[3].y - verts[0].y) / 5};
    Point tp{verts[0].x + dist(mt) * (verts[3].x - verts[0].x), verts[0].y + dist(mt) * (verts[3].y - verts[0].y)};
    double expected = f(tp);
    double actual = 0;
    Point pact;
    for (auto i = 0; i < total; ++i)
        actual += sol[i] * GetShapeFunction(i, tp);
    cout << node <<  "\t" << actual << "\t" << expected << endl;
    auto func = [=](const Point& p)
    {
        return Point(4 * p.x * p.x * p.x, 4 * p.y * p.y * p.y);
    };
    Point exps = func(tp);
    Point act{0, 0};
    for (auto i = 0; i < total; ++i)
        act += sol[i] * GetGradShapeFunction(i, tp);
    cout << node <<  "\t" << act.x << "\t" << exps.x << endl;
    cout << node <<  "\t" << act.y << "\t" << exps.y << endl;*/
    return sol[node];
}


