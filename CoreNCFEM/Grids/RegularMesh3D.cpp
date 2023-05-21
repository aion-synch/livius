#include <stdio.h>
#include "RegularMesh3D.h"
#include "../FiniteElements/Cube.h"
#include "../FiniteElements/Edge.h"
#include <iostream>
#include <algorithm>
#include <numeric>
using namespace std;
using namespace corenc;
using namespace Mesh;

template<class T>
vector<size_t> sort_indexes(const vector<T> &v)
{
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}

CRegularMesh3D::CRegularMesh3D()
{
    m_order = 1;
}
const double CRegularMesh3D::CompSquare(const Point& p1, const Point& p2, const Point& p3) const
{
    return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
}
CRegularMesh3D::~CRegularMesh3D()
{
    if (m_elems.size() > 0)
    {
        auto sz(m_elems.size());
        for (auto i = 0; i < sz; ++i)
            delete m_elems[i];
        vector<CElement<>*>().swap(m_elems);
    }
    if (m_edges.size() > 0)
    {
        auto sz(m_edges.size());
        for (auto i = 0; i < sz; ++i)
            delete m_edges[i];
        vector<CElement<>*>().swap(m_edges);
    }
    if (m_points.size() > 0)
    {
        vector<Point>().swap(m_points);
    }
    if (m_basepoints.size() > 0)
        vector<Point>().swap(m_basepoints);
    if (m_elemsbase.size() > 0)
    {
        auto sz(m_elemsbase.size());
        for (auto i = 0; i < sz; ++i)
            delete m_elemsbase[i];
        vector<CElement<>*>().swap(m_elemsbase);
    }
    if (m_edgesbase.size() > 0)
    {
        auto sz(m_edgesbase.size());
        for (auto i = 0; i < sz; ++i)
            delete m_edgesbase[i];
        vector<CElement<>*>().swap(m_edgesbase);
    }
    if (m_params.size() > 0)
        std::map<int, CParameter>().swap(m_params);
}

CRegularMesh3D::CRegularMesh3D(const string& file_name)
{
    //cout << file_name << endl;
    //cout << "Reading the grid from the file..." << endl;
    ifstream in{ file_name };
    if (!in.is_open())
    {
        //cout << "Error! Couldn't load the grid\n";
        return;
    }
    CFiniteElement<CCube, CCubeBasis, int> ele{};
    CFiniteElement<CEdge, CEdgeLinearBasis, int> edge{};
    CEdge edg;
    CEdgeLinearBasis edglb;
    float temp;
    int num, tp;
    int nums;
    vector<int> nods(3), nds(3);
    unsigned int k{ 0 };
    Point points[3];
    string tmp;
    Point p;
    //EdgeQuad edge;
    //vector<TriangleQuad> trs;
    vector<int> tps;
    auto scal = [&](const Point& p1, const Point& p2) {return p1.x*p2.x + p1.y*p2.y + p1.z * p2.z; };
    try
    {
        in >> tmp >> temp >> temp >> temp >> tmp >> tmp;
        in >> num;
        if (num < 1)
        {
            //cout << "\nError: size is less than 1" << endl;
            throw;
        }
        m_points.resize(num);
        for (int i = 0; i < num; ++i)
            in >> tmp >> m_points[i].x >> m_points[i].y >> tmp;
        in >> tmp >> tmp;
        in >> nums;
        vector<pair<double, double>> vec(num);
        vector<int> idx(num);
        for (int i = 0; i < num; ++i)
        {
            vec[i].second = m_points[i].x;
            vec[i].first = m_points[i].y;
        }
        int iz = 0;
        for (auto i : sort_indexes(vec))
        {
            idx[i] = iz;
            ++iz;
        }
        std::sort(std::begin(vec), std::end(vec));
        for (int i = 0; i < num; ++i)
        {
            m_points[i].x = vec[i].second;
            m_points[i].y = vec[i].first;
        }
        for (int i = 0; i < nums; ++i)
        {
            in >> num >> tp;
            switch (tp)
            {
            case 15:
                in >> tp >> tp >> tp >> tp; break;
            case 1:
                in >> tp >> tp >> tp >> nods[0] >> nods[1];
                --nods[0]; --nods[1];
                nods[0] = idx[nods[0]];
                nods[1] = idx[nods[1]];
                points[0] = m_points[nods[0]];
                points[1] = m_points[nods[1]];
                m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>(&nods[0], points, 2, tp));
                break;
            case 2:
            {
                in >> tp >> tp >> tp >> nods[0] >> nods[1] >> nods[2];
                --nods[0]; --nods[1]; --nods[2];
                nods[0] = idx[nods[0]];
                nods[1] = idx[nods[1]];
                nods[2] = idx[nods[2]];
                std::sort(nods.begin(), nods.end());
                nds[0] = nods[2];
                nds[1] = nods[0];
                points[0] = m_points[nds[0]];
                points[1] = m_points[nds[1]];
                points[2] = m_points[nods[1]];
                k = 0;
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{&nds[0], points, 2, -1} };
                for (; k < m_edges.size(); ++k)
                {
                    if (edg1 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[k]))
                    {
                        if (m_edges[k]->GetNeighbour(0) > -1)
                            m_edges[k]->SetNeighbour(1, (int)m_elems.size());
                        else
                            m_edges[k]->SetNeighbour(0, (int)m_elems.size());
                        break;
                    }
                }
                if (k == m_edges.size())
                {
                    m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg1 });
                    m_edges[m_edges.size() - 1]->SetNeighbour(0, m_elems.size());
                    if (scal(m_edges[m_edges.size() - 1]->GetNormal(), points[2]) > 0)
                        m_edges[m_edges.size() - 1]->ReverseNormal();
                }
                nds[0] = nods[1];
                nds[1] = nods[2];
                points[0] = m_points[nds[0]];
                points[1] = m_points[nds[1]];
                points[2] = m_points[nods[0]];
                k = 0;
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{&nds[0], points, 2, -1} };
                for (; k < m_edges.size(); ++k)
                {
                    if (edg2 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[k]))
                    {
                        if (m_edges[k]->GetNeighbour(0) > -1)
                            m_edges[k]->SetNeighbour(1, (int)m_elems.size());
                        else
                            m_edges[k]->SetNeighbour(0, (int)m_elems.size());
                        break;
                    }
                }
                if (k == m_edges.size())
                {
                    m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg2 });
                    m_edges[m_edges.size() - 1]->SetNeighbour(0, m_elems.size());
                    if (scal(m_edges[m_edges.size() - 1]->GetNormal(), points[2]) > 0)
                        m_edges[m_edges.size() - 1]->ReverseNormal();
                }

                points[0] = m_points[nods[0]];
                points[1] = m_points[nods[1]];
                points[2] = m_points[nods[2]];
                k = 0;
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg3{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{&nods[0], points, 2, -1} };

                for (; k < m_edges.size(); ++k)
                {
                    if (edg3 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[k]))
                    {
                        if (m_edges[k]->GetNeighbour(0) > -1)
                            m_edges[k]->SetNeighbour(1, (int)m_elems.size());
                        else
                            m_edges[k]->SetNeighbour(0, (int)m_elems.size());
                        break;
                    }
                }
                if (k == m_edges.size())
                {
                    m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg3 });
                    m_edges[m_edges.size() - 1]->SetNeighbour(0, m_elems.size());
                    if (scal(m_edges[m_edges.size() - 1]->GetNormal(), points[2]) > 0)
                        m_edges[m_edges.size() - 1]->ReverseNormal();
                }
                std::sort(nods.begin(), nods.end());
                /*if (scal(Point(points[0], points[], m_points[nods[1]]) >= 0)
                {

                }*/
                //points[0] = m_points[nods[0]];
                //points[1] = m_points[nods[1]];
                //points[2] = m_points[nods[2]];
                CFiniteElement<CCube, CCubeBasis> tr{ CFiniteElement<CCube, CCubeBasis>{&nods[0], points, 3, tp} };
                if (scal(tr.GetNormal(), Point(0, 0, 0.5)) < 0)
                {
                    nods[1] = tr.GetNode(2);
                    nods[2] = tr.GetNode(1);
                    points[1] = m_points[nods[1]];
                    points[2] = m_points[nods[2]];
                }
                m_elems.push_back(new CFiniteElement<CCube, CCubeBasis>{ &nods[0], points, 3, tp });

                break;
            }
            /*case 2:
            cin >> tp >> tp >> tp >> nods[0] >> nods[1] >> nods[2];
            --nods[0]; --nods[1]; --nods[2];
            points[0] = m_points[nods[0]];
            points[1] = m_points[nods[1]];
            points[2] = m_points[nods[2]];
            edge.m_nodes[0] = -1;
            edge.m_nodes[1] = nods[0];
            edge.m_nodes[2] = nods[1];
            edge.CompMes(m_points[nods[0]], m_points[nods[1]]);
            edge.m_normal.x = points[0].y - points[1].y;
            edge.m_normal.y = points[1].x - points[0].x;
            //edge.m_normal.x = m_points[nods[0]].z - m_points[nods[1]].z;
            edge.order = 2;
            edge.nfem[0] = trs.size();
            k = 0;
            for (; k < m_edges.size(); ++k)
            {
            if (m_edges[k] == edge)
            {
            m_edges[k].m_normal.x = points[0].y - points[1].y;
            m_edges[k].m_normal.y = points[1].x - points[0].x;
            if (m_edges[k].nfem[0] > -1)
            m_edges[k].nfem[1] = trs.size();
            else
            m_edges[k].nfem[0] = trs.size();
            break;
            }
            }
            if (k == m_edges.size())
            {
            m_edges.push_back(edge);
            if (scal(m_edges[m_edges.size() - 1].GetNormal(), points[2]) > 0)
            {
            m_edges[m_edges.size() - 1].m_normal.x = -m_edges[m_edges.size() - 1].m_normal.x;
            m_edges[m_edges.size() - 1].m_normal.y = -m_edges[m_edges.size() - 1].m_normal.y;
            }
            }
            trd.n1 = k;
            edge.m_nodes[0] = -1;
            edge.m_nodes[1] = nods[1];
            edge.m_nodes[2] = nods[2];
            edge.CompMes(m_points[nods[1]], m_points[nods[2]]);
            edge.m_normal.x = points[1].y - points[2].y;
            edge.m_normal.y = points[2].x - points[1].x;
            edge.nfem[0] = trs.size();
            k = 0;
            for (; k < m_edges.size(); ++k)
            {
            if (m_edges[k] == edge)
            {
            m_edges[k].m_normal.x = points[1].y - points[2].y;
            m_edges[k].m_normal.y = points[2].x - points[1].x;
            if (m_edges[k].nfem[0] > -1)
            m_edges[k].nfem[1] = trs.size();
            else
            m_edges[k].nfem[0] = trs.size();
            break;
            }
            }
            if (k == m_edges.size())
            {
            m_edges.push_back(edge);
            if (scal(m_edges[m_edges.size() - 1].GetNormal(), points[0]) > 0)
            {
            m_edges[m_edges.size() - 1].m_normal.x = -m_edges[m_edges.size() - 1].m_normal.x;
            m_edges[m_edges.size() - 1].m_normal.y = -m_edges[m_edges.size() - 1].m_normal.y;
            }
            }
            trd.n2 = k;
            edge.m_nodes[0] = -1;
            edge.m_nodes[1] = nods[2];
            edge.m_nodes[2] = nods[0];
            edge.CompMes(m_points[nods[2]], m_points[nods[0]]);
            edge.m_normal.x = points[2].y - points[0].y;
            edge.m_normal.y = points[0].x - points[2].x;
            edge.nfem[0] = trs.size();
            k = 0;
            for (; k < m_edges.size(); ++k)
            {
            if (m_edges[k] == edge)
            {
            m_edges[k].m_normal.x = points[2].y - points[0].y;
            m_edges[k].m_normal.y = points[0].x - points[2].x;
            if (m_edges[k].nfem[0] > -1)
            m_edges[k].nfem[1] = trs.size();
            else
            m_edges[k].nfem[0] = trs.size();
            break;
            }
            }
            if (k == m_edges.size())
            {
            m_edges.push_back(edge);
            if (scal(m_edges[m_edges.size() - 1].GetNormal(), points[1]) > 0)
            {
            m_edges[m_edges.size() - 1].m_normal.x = -m_edges[m_edges.size() - 1].m_normal.x;
            m_edges[m_edges.size() - 1].m_normal.y = -m_edges[m_edges.size() - 1].m_normal.y;
            }
            }
            trd.n3 = k;
            //if (get<6>(m_params.find(tp)->second) == 2)
            //m_fems.push_back(new CFiniteElement < TriangleQuad, QuadBasisTr >{ nods, points, tp });
            //else
            //m_fems.push_back(new CFiniteElement < TriangleQuad, QuadBasisTr >{ nods, points, tp });
            trds.push_back(trd);
            trs.push_back(TriangleQuad{ nods });
            tps.push_back(tp);
            //m_fems.push_back(new CFiniteElement < Triangle, LinearBasisTr > {nods, points, tp});
            break;*/
            default:
                break;
            }
        }
    }
    catch (exception &e)
    {
        //cout << "Sorry, error: ";
        //cout << e.what() << endl;
        exit(0);
    }
    m_basepoints = m_points;
    const auto& sz_el = m_elems.size();
    const auto& sz_ed = m_edges.size();
    m_elemsbase.resize(sz_el);
    m_edgesbase.resize(sz_ed);
    for (int i = 0; i < sz_el; ++i)
    {
        nods[0] = m_elems[i]->GetNode(0);
        nods[1] = m_elems[i]->GetNode(1);
        nods[2] = m_elems[i]->GetNode(2);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        m_elemsbase[i] = new CFiniteElement<CCube, CCubeBasis>{ &nods[0], points, 3, m_elems[i]->GetType() };
    }
    for (int i = 0; i < sz_ed; ++i)
    {
        nods[0] = m_edges[i]->GetNode(0);
        nods[1] = m_edges[i]->GetNode(1);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        m_edgesbase[i] = new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ &nods[0], points, 2, m_edges[i]->GetType() };
    }
    m_order = 1;
}

CRegularMesh3D::CRegularMesh3D(const CRegularMesh3D &tr)
{
    const int sz_el = tr.m_elems.size();
    const int sz_pt = tr.m_points.size();
    const int sz_bpt = tr.m_basepoints.size();
    const int sz_bel = tr.m_elemsbase.size();
    const int sz_ed = tr.m_edges.size();
    const int sz_bed = tr.m_edgesbase.size();
    m_elems.resize(sz_el);
    m_edges.resize(sz_ed);
    m_points.resize(sz_pt);
    m_basepoints.resize(sz_bpt);
    m_elemsbase.resize(sz_bel);
    m_edgesbase.resize(sz_bed);
    int i = 0;
    for (i = 0; i < sz_el; ++i)
        m_elems[i] = tr.m_elems[i]->Clone();
    for (i = 0; i < sz_ed; ++i)
        m_edges[i] = tr.m_edges[i]->Clone();
    for (i = 0; i < sz_pt; ++i)
        m_points[i] = tr.m_points[i];
    for (i = 0; i < sz_bpt; ++i)
        m_basepoints[i] = tr.m_basepoints[i];
    for (i = 0; i < sz_bel; ++i)
        m_elemsbase[i] = tr.m_elemsbase[i]->Clone();
    for (i = 0; i < sz_bed; ++i)
        m_edgesbase[i] = tr.m_edgesbase[i]->Clone();
    m_params = tr.m_params;
    m_order = tr.m_order;
    m_bnds.resize(tr.m_bnds.size());
    for (i = 0; i < m_bnds.size(); ++i)
        m_bnds[i] = tr.m_bnds[i];
    m_offsets = tr.m_offsets;
    m_inodes = tr.m_inodes;
}

CRegularMesh3D::CRegularMesh3D(const double x1, const double y1, const double x2, const double y2, const int nx, const int ny)
{
    const double hx = (x2 - x1) / nx;
    const double hy = (y2 - y1) / ny;
    const size_t size_p = (nx + 1) * (ny + 1);
    const size_t size_f = nx * ny;
    m_points.resize(size_p);
    for (size_t i = 0; i < ny + 1; ++i)
        for (size_t j = 0; j < nx + 1; ++j)
            m_points[i * (nx + 1) + j] = Point(x1 + hx * j, y1 + hy * i);
    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            Point temp[4];
            int nodes[4];
            nodes[0] = i * (nx + 1) + j;
            nodes[1] = i * (nx + 1) + j + 1;
            nodes[2] = (i + 1) * (nx + 1) + j;
            nodes[3] = (i + 1) * (nx + 1) + j + 1;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            temp[2] = m_points[nodes[2]];
            temp[3] = m_points[nodes[3]];

            m_elems.push_back(new CFiniteElement<CCube, CCubeBasis>{ nodes, temp, 4, 0 });
        }
    }
    for (size_t i = 0; i < nx; ++i)
    {
        Point temp[2];
        int nodes[2];
        temp[0] = m_points[i];
        temp[1] = m_points[i + 1];
        nodes[0] = i;
        nodes[1] = i + 1;
        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 0 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, i);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        m_edges[m_edges.size() - 1]->ReverseNormal();
        nodes[0] = ny * (nx + 1) + i;
        nodes[1] = ny * (nx + 1) + i + 1;
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        //m_edges[m_edges.size() - 1]->ReverseNormal();

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 1 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, (ny - 1) * nx + i);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
    }

    for (size_t i = 0; i < ny; ++i)
    {
        Point temp[2];
        int nodes[2];
        nodes[0] = i * (nx + 1);
        nodes[1] = (i + 1) * (nx + 1);
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 2 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, i * nx);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        m_edges[m_edges.size() - 1]->ReverseNormal();

        nodes[0] = nx + i * (nx + 1);
        nodes[1] = nx + (i + 1) * (nx + 1);
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        //m_edges[m_edges.size() - 1]->ReverseNormal();

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 3 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, nx - 1 + i * nx);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, - 1);
    }
    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx - 1; ++j)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = 1 + j + (nx + 1) * i;
            nodes[1] = 1 + j + (nx + 1) * (i + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, j + i * nx + 1);
        }
    }

    for (size_t i = 0; i < ny - 1; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = nx + 1 + j + (nx + 1) * i;
            nodes[1] = nx + 2 + j + (nx + 1) * i;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, j + (i + 1) * nx);
        }
    }
}

CRegularMesh3D::CRegularMesh3D(const Point& p1, const Point& p2, const int nx, const int ny)
{
    const auto x1 = p1.x;
    const auto y1 = p1.y;
    const auto x2 = p2.x;
    const auto y2 = p2.y;
    const double hx = (x2 - x1) / nx;
    const double hy = (y2 - y1) / ny;
    const size_t size_p = (nx + 1) * (ny + 1);
    const size_t size_f = nx * ny;
    m_offsets.resize(size_p);
    m_bnds.resize(size_p);
    m_points.resize(size_p);
    m_inodes = 2 * nx + 2 * ny;
    for (int k = 0; k < size_p; ++k)
        m_offsets[k] -= nx + 2;
    for (size_t i = 0; i < ny + 1; ++i)
        for (size_t j = 0; j < nx + 1; ++j)
            m_points[i * (nx + 1) + j] = Point(x1 + hx * j, y1 + hy * i);
    for (int i = 1; i < ny; ++i)
    {
        for (int k = (i + 1) * (nx + 1); k < size_p; ++k)
            m_offsets[k] -= 2;
    }
    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            Point temp[4];
            int nodes[4];
            nodes[0] = i * (nx + 1) + j;
            nodes[1] = i * (nx + 1) + j + 1;
            nodes[2] = (i + 1) * (nx + 1) + j;
            nodes[3] = (i + 1) * (nx + 1) + j + 1;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            temp[2] = m_points[nodes[2]];
            temp[3] = m_points[nodes[3]];

            m_elems.push_back(new CFiniteElement<CCube, CCubeBasis>{ nodes, temp, 4, 0 });
        }
    }
    for (size_t i = 0; i < nx; ++i)
    {
        Point temp[2];
        int nodes[2];
        temp[0] = m_points[i];
        temp[1] = m_points[i + 1];
        nodes[0] = i;
        nodes[1] = i + 1;
        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;
        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 0 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, i);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        m_edges[m_edges.size() - 1]->ReverseNormal();
        nodes[0] = ny * (nx + 1) + i;
        nodes[1] = ny * (nx + 1) + i + 1;
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        //m_edges[m_edges.size() - 1]->ReverseNormal();

        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;
        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 1 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, (ny - 1) * nx + i);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
    }

    for (size_t i = 0; i < ny; ++i)
    {
        Point temp[2];
        int nodes[2];
        nodes[0] = i * (nx + 1);
        nodes[1] = (i + 1) * (nx + 1);
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 2 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, i * nx);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        m_edges[m_edges.size() - 1]->ReverseNormal();

        nodes[0] = nx + i * (nx + 1);
        nodes[1] = nx + (i + 1) * (nx + 1);
        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        //m_edges[m_edges.size() - 1]->ReverseNormal();

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 3 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, nx - 1 + i * nx);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, - 1);
    }
    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx - 1; ++j)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = 1 + j + (nx + 1) * i;
            nodes[1] = 1 + j + (nx + 1) * (i + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, j + i * nx + 1);
        }
    }

    for (size_t i = 0; i < ny - 1; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = nx + 1 + j + (nx + 1) * i;
            nodes[1] = nx + 2 + j + (nx + 1) * i;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, j + (i + 1) * nx);
        }
    }
}

CRegularMesh3D::CRegularMesh3D(const Point& p1, const Point& p2, const int nx, const int ny, const int px, const int py)
{
    /*if ((px == 2) && (py == 1))
    {
        const auto x1 = p1.x;
        const auto y1 = p1.y;
        const auto x2 = p2.x;
        const auto y2 = p2.y;
        const double hx = (x2 - x1) / nx;
        const double hy = (y2 - y1) / ny;
        const size_t size_p = (nx + 1) * (ny + 1);
        const size_t size_f = nx * ny;
        const size_t size_p2 = nx;
        m_points.resize(size_p + nx * (ny + 1));
        for (size_t i = 0; i < ny + 1; ++i)
            for (size_t j = 0; j < nx + 1; ++j)
                m_points[i * (nx + 1) + j] = Point(x1 + hx * j, y1 + hy * i);
        for (auto i = 0; i < ny + 1; ++i)
            for (auto j = 0; j < nx; ++j)
                m_points[(ny + 1) * (nx + 1) + i * nx + j] = Point(x1 + hx / 2 + hx * j, y1 + hy * i);
        for (size_t i = 0; i < ny; ++i)
        {
            for (size_t j = 0; j < nx; ++j)
            {
                Point temp[6];
                int nodes[6];
                nodes[0] = i * (nx + 1) + j;
                nodes[1] = i * (nx + 1) + j + 1;
                nodes[2] = (i + 1) * (nx + 1) + j;
                nodes[3] = (i + 1) * (nx + 1) + j + 1;
                nodes[4] = (ny + 1) * (nx + 1) + i * (nx) + j;
                nodes[5] = (ny + 1) * (nx + 1) + (i + 1) * (nx) + j;
                temp[0] = m_points[nodes[0]];
                temp[1] = m_points[nodes[1]];
                temp[2] = m_points[nodes[2]];
                temp[3] = m_points[nodes[3]];
                temp[4] = m_points[nodes[4]];
                temp[5] = m_points[nodes[5]];

                m_elems.push_back(new CFiniteElement<CCube, CCubeBasis2x>{ nodes, temp, 6, 0 });
                //for (int ii = 0; ii < m_elems[m_elems.size() - 1]->GetDoFs(); ++ii)
                  //  cout << m_elems[m_elems.size() - 1]->GetNode(ii) << "\t";
                //cout << endl;
            }
        }
        for (size_t i = 0; i < nx; ++i)
        {
            Point temp[2];
            int nodes[2];
            temp[0] = m_points[i];
            temp[1] = m_points[i + 1];
            nodes[0] = i;
            nodes[1] = i + 1;
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 0 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, i);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
            m_edges[m_edges.size() - 1]->ReverseNormal();
            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(2, (ny + 1) * (nx + 1) + i);
            nodes[0] = ny * (nx + 1) + i;
            nodes[1] = ny * (nx + 1) + i + 1;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            //m_edges[m_edges.size() - 1]->ReverseNormal();

            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, (ny - 1) * nx + i);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(2, (ny + 1) * (nx + 1) + ny * (nx) + i);
        }

        for (size_t i = 0; i < ny; ++i)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = i * (nx + 1);
            nodes[1] = (i + 1) * (nx + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];

            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 2 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
            m_edges[m_edges.size() - 1]->ReverseNormal();

            nodes[0] = nx + i * (nx + 1);
            nodes[1] = nx + (i + 1) * (nx + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            //m_edges[m_edges.size() - 1]->ReverseNormal();

            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 3 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, nx - 1 + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, - 1);
        }

        // y axis | edges
        for (size_t i = 0; i < ny; ++i)
        {
            for (size_t j = 0; j < nx - 1; ++j)
            {
                Point temp[2];
                int nodes[2];
                nodes[0] = 1 + j + (nx + 1) * i;
                nodes[1] = 1 + j + (nx + 1) * (i + 1);
                temp[0] = m_points[nodes[0]];
                temp[1] = m_points[nodes[1]];
                m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
                m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
                m_edges[m_edges.size() - 1]->SetNeighbour(1, j + i * nx + 1);
            }
        }

        // x axis - edges
        for (size_t i = 0; i < ny - 1; ++i)
        {
            for (size_t j = 0; j < nx; ++j)
            {
                Point temp[2];
                int nodes[2];
                nodes[0] = nx + 1 + j + (nx + 1) * i;
                nodes[1] = nx + 2 + j + (nx + 1) * i;
                temp[0] = m_points[nodes[0]];
                temp[1] = m_points[nodes[1]];
                m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
                m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
                m_edges[m_edges.size() - 1]->SetNeighbour(1, j + (i + 1) * nx);
                //m_edges[m_edges.size() - 1]->IncreaseOrder();
                //m_edges[m_edges.size() - 1]->SetNode(2, (ny + 1) * (nx + 1) + (i + 1) * nx + j);
            }
        }
    }
    else if ((py == 2) && (px == 1))
    {
        const auto x1 = p1.x;
        const auto y1 = p1.y;
        const auto x2 = p2.x;
        const auto y2 = p2.y;
        const double hx = (x2 - x1) / nx;
        const double hy = (y2 - y1) / ny;
        const size_t size_p = (nx + 1) * (ny + 1);
        const size_t size_f = nx * ny;
        const size_t size_p2 = nx;
        m_points.resize(size_p + (nx + 1) * ny);
        for (size_t i = 0; i < ny + 1; ++i)
            for (size_t j = 0; j < nx + 1; ++j)
                m_points[i * (nx + 1) + j] = Point(x1 + hx * j, y1 + hy * i);
        for (auto i = 0; i < nx + 1; ++i)
            for (auto j = 0; j < ny; ++j)
                m_points[(ny + 1) * (nx + 1) + i * ny + j] = Point(x1 + hx * i, y1 + hy / 2 + hy * j);
        cout << "pts\n";
        for (auto & it : m_points)
            cout << it.x << "\t" << it.y << endl;
        for (size_t i = 0; i < ny; ++i)
        {
            for (size_t j = 0; j < nx; ++j)
            {
                Point temp[6];
                int nodes[6];
                nodes[0] = i * (nx + 1) + j;
                nodes[1] = i * (nx + 1) + j + 1;
                nodes[2] = (i + 1) * (nx + 1) + j;
                nodes[3] = (i + 1) * (nx + 1) + j + 1;
                nodes[4] = (ny + 1) * (nx + 1) + j * (ny) + i;
                nodes[5] = (ny + 1) * (nx + 1) + (j + 1) * (ny) + i;
                temp[0] = m_points[nodes[0]];
                temp[1] = m_points[nodes[1]];
                temp[2] = m_points[nodes[2]];
                temp[3] = m_points[nodes[3]];
                temp[4] = m_points[nodes[4]];
                temp[5] = m_points[nodes[5]];

                m_elems.push_back(new CFiniteElement<CCube, CCubeBasis2y>{ nodes, temp, 6, 0 });
                //for (int ii = 0; ii < m_elems[m_elems.size() - 1]->GetDoFs(); ++ii)
                  //  cout << m_elems[m_elems.size() - 1]->GetNode(ii) << "\t";
                //cout << endl;
            }
        }
        for (size_t i = 0; i < nx; ++i)
        {
            Point temp[2];
            int nodes[2];
            temp[0] = m_points[i];
            temp[1] = m_points[i + 1];
            nodes[0] = i;
            nodes[1] = i + 1;
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 0 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, i);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
            m_edges[m_edges.size() - 1]->ReverseNormal();
            nodes[0] = ny * (nx + 1) + i;
            nodes[1] = ny * (nx + 1) + i + 1;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            //m_edges[m_edges.size() - 1]->ReverseNormal();

            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, (ny - 1) * nx + i);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        }

        for (size_t i = 0; i < ny; ++i)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = i * (nx + 1);
            nodes[1] = (i + 1) * (nx + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];

            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 2 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
            m_edges[m_edges.size() - 1]->ReverseNormal();
            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(2, (ny + 1) * (nx + 1) + i);

            nodes[0] = nx + i * (nx + 1);
            nodes[1] = nx + (i + 1) * (nx + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            //m_edges[m_edges.size() - 1]->ReverseNormal();

            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 3 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, nx - 1 + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, - 1);
            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(2, (ny + 1) * (nx + 1) + nx * ny + i);
        }

        // y axis | edges
        for (size_t i = 0; i < ny; ++i)
        {
            for (size_t j = 0; j < nx - 1; ++j)
            {
                Point temp[2];
                int nodes[2];
                nodes[0] = 1 + j + (nx + 1) * i;
                nodes[1] = 1 + j + (nx + 1) * (i + 1);
                temp[0] = m_points[nodes[0]];
                temp[1] = m_points[nodes[1]];
                m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
                m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
                m_edges[m_edges.size() - 1]->SetNeighbour(1, j + i * nx + 1);
            }
        }

        // x axis - edges
        for (size_t i = 0; i < ny - 1; ++i)
        {
            for (size_t j = 0; j < nx; ++j)
            {
                Point temp[2];
                int nodes[2];
                nodes[0] = nx + 1 + j + (nx + 1) * i;
                nodes[1] = nx + 2 + j + (nx + 1) * i;
                temp[0] = m_points[nodes[0]];
                temp[1] = m_points[nodes[1]];
                m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
                m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
                m_edges[m_edges.size() - 1]->SetNeighbour(1, j + (i + 1) * nx);
                m_edges[m_edges.size() - 1]->IncreaseOrder();
                m_edges[m_edges.size() - 1]->SetNode(2, (ny + 1) * (nx + 1) + (i + 1) * nx + j);
            }
        }
    }
    else
    {*/
    const auto x1 = p1.x;
    const auto y1 = p1.y;
    const auto x2 = p2.x;
    const auto y2 = p2.y;
    const double hx = (x2 - x1) / nx;
    const double hy = (y2 - y1) / ny;
    //const double hpx = hx / px;
    //const double hpy = hy / py;
    const size_t size_p = (nx + 1) * (ny + 1);
    //int size_p = (px * nx + 1) * (py * ny + 1);
    const size_t size_f = nx * ny;
    int total_degree = (px + 1) * (py + 1);
    m_points.resize(size_p);
    int npy = ny * py + 1;
    int npx = nx * px + 1;
    for (size_t i = 0; i < ny + 1; ++i)
        for (size_t j = 0; j < nx + 1; ++j)
            m_points[i * (nx + 1) + j] = Point(x1 + hx * j, y1 + hy * i);
    /*for (int ii = 0; ii < py; ++ii)
        for (int jj = 0; jj < px; ++jj)
        {
            for (size_t i = 0; i < ny + 1; ++i)
                for (size_t j = 0; j < nx + 1; ++j)
                {
                    m_points[ii * (px * py * nx * ny + 1) + jj * (nx + 1) * (ny + 1) + i * (nx + 1) + j] = Point(x1 + hx * j / (jj + 1), y1 + hy * i / (ii + 1));
                    cout << ii << "\t" << jj << "\t" << i << "\t" << j << endl;
                    cout << ii * (py * ny * px * nx + 1) + jj * (nx + 1) * (ny + 1) + i * (nx + 1) + j << endl;
                }
            cout << endl;
        }*/
    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            Point temp[4];
            int nodes[4];
            nodes[0] = i * (nx + 1) + j;
            nodes[1] = i * (nx + 1) + j + 1;
            nodes[2] = (i + 1) * (nx + 1) + j;
            nodes[3] = (i + 1) * (nx + 1) + j + 1;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            temp[2] = m_points[nodes[2]];
            temp[3] = m_points[nodes[3]];
            //cout << "cr" << endl;
            m_elems.push_back(new CFiniteElement<CCube, CCubeBasis>{ nodes, temp, 4, 0 });

            //m_elems[m_elems.size() - 1]->SetOrder(px, py);
            //cout << m_elems[m_elems.size() - 1]->GetDoFs() << endl;
            for (int ii = 1; ii < py; ++ii)
            {
                for (int jj = 1; jj < px; ++jj)
                {
                    m_points.push_back(Point(temp[0].x + (temp[3].x - temp[0].x) / (jj + 1), temp[0].y + (temp[3].y - temp[0].y) / (ii + 1)));
                    m_elems[m_elems.size() - 1]->SetNode(3 + 2 * (px - 1) + 2 * (py - 1) + (ii - 1) * (px - 1) + jj, m_points.size() - 1);
                }
            }
        }
    }
    for (size_t i = 0; i < nx; ++i)
    {
        Point temp[2];
        int nodes[2];
        temp[0] = m_points[i];
        temp[1] = m_points[i + 1];
        nodes[0] = i;
        nodes[1] = i + 1;
        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 0 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, i);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        m_edges[m_edges.size() - 1]->ReverseNormal();

        for (int j = 1; j < px; ++j)
        {
            m_points.push_back(Point(temp[0].x + (temp[1].x - temp[0].x) / (j + 1), temp[0].y));
            m_elems[i]->SetNode(3 + j, m_points.size() - 1);
            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(1 + j, m_points.size() - 1);
        }

        nodes[0] = ny * (nx + 1) + i;
        nodes[1] = ny * (nx + 1) + i + 1;
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        //m_edges[m_edges.size() - 1]->ReverseNormal();

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 1 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, (ny - 1) * nx + i);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        for (int j = 1; j < px; ++j)
        {
            m_points.push_back(Point(temp[0].x + (temp[1].x - temp[0].x) / (j + 1), temp[0].y));
            m_elems[(ny - 1) * nx + i]->SetNode(3 + px + j - 1, m_points.size() - 1);

            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(1 + j, m_points.size() - 1);
        }
    }

    for (size_t i = 0; i < ny; ++i)
    {
        Point temp[2];
        int nodes[2];
        nodes[0] = i * (nx + 1);
        nodes[1] = (i + 1) * (nx + 1);
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 2 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, i * nx);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
        m_edges[m_edges.size() - 1]->ReverseNormal();

        for (int j = 1; j < py; ++j)
        {
            m_points.push_back(Point(temp[0].x, temp[0].y + (temp[1].y - temp[0].y) / (j + 1)));
            m_elems[i * nx]->SetNode(3 + 2 * (px - 1) + j, m_points.size() - 1);

            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(1 + j, m_points.size() - 1);
        }

        nodes[0] = nx + i * (nx + 1);
        nodes[1] = nx + (i + 1) * (nx + 1);
        temp[0] = m_points[nodes[0]];
        temp[1] = m_points[nodes[1]];
        //m_edges[m_edges.size() - 1]->ReverseNormal();

        m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 3 });
        m_edges[m_edges.size() - 1]->SetNeighbour(0, nx - 1 + i * nx);
        m_edges[m_edges.size() - 1]->SetNeighbour(1, - 1);
        for (int j = 1; j < py; ++j)
        {
            m_points.push_back(Point(temp[0].x, temp[0].y + (temp[1].y - temp[0].y) / (j + 1)));
            m_elems[nx - 1 + i * nx]->SetNode(3 + 2 * (px - 1) + py - 1 + j, m_points.size() - 1);

            m_edges[m_edges.size() - 1]->IncreaseOrder();
            m_edges[m_edges.size() - 1]->SetNode(1 + j, m_points.size() - 1);
        }
    }
    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx - 1; ++j)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = 1 + j + (nx + 1) * i;
            nodes[1] = 1 + j + (nx + 1) * (i + 1);
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, j + i * nx + 1);

            for (int jj = 1; jj < py; ++jj)
            {
                m_points.push_back(Point(temp[0].x, temp[0].y + (temp[1].y - temp[0].y) / (jj + 1)));
                m_elems[j + i * nx]->SetNode(3 + 2 * (px - 1) + py + jj - 1, m_points.size() - 1);
                m_elems[j + i * nx + 1]->SetNode(3 + 2 * (px - 1) + jj, m_points.size() - 1);

                m_edges[m_edges.size() - 1]->IncreaseOrder();
                m_edges[m_edges.size() - 1]->SetNode(1 + jj, m_points.size() - 1);
            }
        }
    }

    for (size_t i = 0; i < ny - 1; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            Point temp[2];
            int nodes[2];
            nodes[0] = nx + 1 + j + (nx + 1) * i;
            nodes[1] = nx + 2 + j + (nx + 1) * i;
            temp[0] = m_points[nodes[0]];
            temp[1] = m_points[nodes[1]];
            m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
            m_edges[m_edges.size() - 1]->SetNeighbour(0, j + i * nx);
            m_edges[m_edges.size() - 1]->SetNeighbour(1, j + (i + 1) * nx);

            for (int jj = 1; jj < px; ++jj)
            {
                m_points.push_back(Point(temp[0].x + (temp[1].x - temp[0].x) / (jj + 1), temp[0].y));
                m_elems[j + i * nx]->SetNode(3 + px + jj - 1, m_points.size() - 1);
                m_elems[j + (i + 1) * nx]->SetNode(3 + jj, m_points.size() - 1);

                m_edges[m_edges.size() - 1]->IncreaseOrder();
                m_edges[m_edges.size() - 1]->SetNode(1 + jj, m_points.size() - 1);
            }
        }
    }
}

const int CRegularMesh3D::refine_h()
{
    return 1;
}

const int CRegularMesh3D::refine_p()
{
    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    return 1;
}

const unsigned int CRegularMesh3D::GetNumberOfNodes() const
{
    return (unsigned int)m_points.size();
}

const Mesh::Point CRegularMesh3D::GetNode(const unsigned int n) const
{
    return m_points[n];
}
const unsigned int CRegularMesh3D::GetNumberOfElements() const
{
    return (unsigned int)m_elems.size();
}

const int CRegularMesh3D::interpolate(const int node) const
{
    /*for (int i = 0; i < m_offsets.size(); ++i)
        cout << m_offsets[i] << endl;*/
    //if (node == 9)
        //cout << 9 << endl;
    if (m_bnds[node])
        return -1;
    return node + m_offsets[node];
}

const int CRegularMesh3D::GetNumberOfINodes() const
{
    //auto n = sqrt(m_points.size());
    //return n * n - 4 * n + 4;
    return m_points.size() - m_inodes;
}

const int CRegularMesh3D::FindElement(const Point& test) const
{
    Point p[3];
    double volElem, mvol;
    for (unsigned int i = 0; i < m_elems.size(); ++i)
    {
        const auto& elem = m_elems[i];
        if(test.x >= m_points[elem->GetNode(0)].x && test.x <= m_points[elem->GetNode(1)].x &&
            test.y >= m_points[elem->GetNode(0)].y && test.y <= m_points[elem->GetNode(3)].y)
            return i;
    }
    return -1;
}

const CElement<>* CRegularMesh3D::GetElement(const unsigned int n) const
{
    if (n < m_elems.size())
        return m_elems[n];
    return nullptr;
}

const CElement<>* CRegularMesh3D::GetBoundary(const unsigned int n) const
{
    if (n < m_edges.size())
        return m_edges[n];
    return nullptr;
}

const double CRegularMesh3D::getSolution(const unsigned int element, const unsigned int node) const
{
    return 0.0;
}

const int CRegularMesh3D::updateSolution(const unsigned int element, const unsigned int node, const double value)
{
    return 0;
}

const std::vector<double> CRegularMesh3D::getSolution() const
{
    return std::vector<double>();
}

const int CRegularMesh3D::updateSolution(const std::vector<double>&)
{
    return 0;
}

const int CRegularMesh3D::updateSolution(const unsigned int element, const unsigned int node, CSolution * value)
{
    return 0;
}

const double CRegularMesh3D::getParameter(Parameters param, const unsigned int l, const Point & p) const
{
    if (l > m_elems.size())
        return 0.;
    const auto& type = m_elems[l]->GetType();
    switch (param)
    {
    case corenc::Parameters::DIFFUSION:
        return m_params.find(type)->second.GetDiffusion(p);
        break;
    case corenc::Parameters::MASS:
        return m_params.find(type)->second.GetMass(p);
        break;
    case corenc::Parameters::ADVECTION:
        return m_params.find(type)->second.GetAdvection(p);
        break;
    default:
        break;
    }
    return 0.;
}

const double CRegularMesh3D::getParameter(Parameters param, const unsigned int l, const int i) const
{
    if (l > m_elems.size())
        return 0.;
    const auto& p = m_points[m_elems[l]->GetNode(i)];
    const auto& type = m_elems[l]->GetType();
    if (m_params.find(type) == m_params.end())
        return 0.;
    switch (param)
    {
    case corenc::Parameters::DIFFUSION:
        return m_params.find(type)->second.GetDiffusion(p);
        break;
    case corenc::Parameters::MASS:
        return m_params.find(type)->second.GetMass(p);
        break;
    case corenc::Parameters::ADVECTION:
        return m_params.find(type)->second.GetAdvection(p);
        break;
    default:
        break;
    }
    return 0.;
}

const int CRegularMesh3D::setParameter(Parameters param, const double, const unsigned int)
{
    return 0;
}

const int CRegularMesh3D::setParameter(const CParameter& p, const unsigned int type)
{
    m_params[type] = p;
    return 0;
}

const int CRegularMesh3D::updateSolution(const unsigned int node, const double value)
{
    return 0;
}


const unsigned int CRegularMesh3D::GetNumberOfBoundaries() const
{
    return (unsigned int)m_edges.size();
}
