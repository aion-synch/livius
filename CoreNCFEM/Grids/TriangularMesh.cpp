#include <stdio.h>
#include "TriangularMesh.h"
#include "../FiniteElements/Triangle.h"
#include "../FiniteElements/Edge.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
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

CTriangularMesh::CTriangularMesh()
{
	m_order = 1;
}
const double CTriangularMesh::CompSquare(const Point& p1, const Point& p2, const Point& p3) const
{
	return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
}
CTriangularMesh::~CTriangularMesh()
{
	if(m_elems.size() > 0)
	{
		auto sz(m_elems.size());
		for(auto i = 0; i < sz; ++i)
			delete m_elems[i];
		vector<CElement<>*>().swap(m_elems);
	}
	if(m_edges.size() > 0)
	{
		auto sz(m_edges.size());
		for(auto i = 0; i < sz; ++i)
			delete m_edges[i];
		vector<CElement<>*>().swap(m_edges);
	}
	if(m_points.size() > 0)
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

CTriangularMesh::CTriangularMesh(const string& file_name)
{
    //cout << file_name << endl;
    //cout << "Reading the grid from the file..." << endl;
    ifstream in{file_name};
    if (!in.is_open())
    {
        //cout << "Error! Couldn't load the grid\n";
        return;
    }
	CFiniteElement<CTriangle, CTriangleBasis, int> ele{};
    CFiniteElement<CTriangle, CTriangleLagrangeBasis, int> ele2{};
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
    auto scal = [&](const Point& p1, const Point& p2){return p1.x*p2.x + p1.y*p2.y + p1.z * p2.z; };
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
			idx[i] = i;
		/*for (int i = 0; i < num; ++i)
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
		}*/
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
					CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{CFiniteElement<CEdge, CEdgeLinearBasis, int>{&nds[0], points, 2, -1}};
					for(; k < m_edges.size(); ++k)
					{
						if(edg1 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[k]))
						{
							if(m_edges[k]->GetNeighbour(0) > -1)
								m_edges[k]->SetNeighbour(1, (int)m_elems.size());
							else
								m_edges[k]->SetNeighbour(0, (int)m_elems.size());
							break;
						}
					}
					if(k == m_edges.size())
					{
						m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{edg1});
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
					CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{CFiniteElement<CEdge, CEdgeLinearBasis, int>{&nds[0], points, 2, -1}};
					for(; k < m_edges.size(); ++k)
					{
						if(edg2 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[k]))
						{
							if(m_edges[k]->GetNeighbour(0) > -1)
								m_edges[k]->SetNeighbour(1, (int)m_elems.size());
							else
								m_edges[k]->SetNeighbour(0, (int)m_elems.size());
							break;
						}
					}
					if(k == m_edges.size())
					{
						m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{edg2});
						m_edges[m_edges.size() - 1]->SetNeighbour(0, m_elems.size());
						if (scal(m_edges[m_edges.size() - 1]->GetNormal(), points[2]) > 0)
							m_edges[m_edges.size() - 1]->ReverseNormal();
					}
					
					points[0] = m_points[nods[0]];
					points[1] = m_points[nods[1]];
					points[2] = m_points[nods[2]];
					k = 0;
					CFiniteElement<CEdge, CEdgeLinearBasis, int> edg3{CFiniteElement<CEdge, CEdgeLinearBasis, int>{&nods[0], points, 2, -1}};
					
					for(; k < m_edges.size(); ++k)
					{
						if(edg3 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[k]))
						{
							if(m_edges[k]->GetNeighbour(0) > -1)
								m_edges[k]->SetNeighbour(1, (int)m_elems.size());
							else
								m_edges[k]->SetNeighbour(0, (int)m_elems.size());
							break;
						}
					}
					if(k == m_edges.size())
					{
						m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{edg3});
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
					//CFiniteElement<CTriangle, CTriangleBasis> tr{ CFiniteElement<CTriangle, CTriangleBasis>{&nods[0], points, 3, tp} };
					//if (scal(tr.GetNormal(), Point(0, 0, 0.5)) < 0)
					//{
					//	nods[1] = tr.GetNode(2);
					//	nods[2] = tr.GetNode(1);
					//	points[1] = m_points[nods[1]];
					//	points[2] = m_points[nods[2]];
					//}
					m_elems.push_back(new CFiniteElement<CTriangle, CTriangleBasis>{&nods[0], points, 1, tp});
					
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
		m_elemsbase[i] = new CFiniteElement<CTriangle, CTriangleBasis>{ &nods[0], points, 1, m_elems[i]->GetType() };
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

CTriangularMesh::CTriangularMesh(const CTriangularMesh &tr)
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
    m_bnds.resize(tr.m_bnds.size());
    for (i = 0; i < m_bnds.size(); ++i)
        m_bnds[i] = tr.m_bnds[i];
    m_offsets = tr.m_offsets;
    m_order = tr.m_order;
}

CTriangularMesh::CTriangularMesh(const Point& p1, const Point& p2, const int nx, const int ny)
{
	const double hx = (p2.x - p1.x) / nx;
	const double hy = (p2.y - p1.y) / ny;
	const size_t size_p = (nx + 1) * (ny + 1);
	const size_t size_f = nx * ny;
	m_points.resize(size_p);
    m_offsets.resize(size_p);
    m_bnds.resize(size_p);
    for (int k = 0; k < size_p; ++k)
        m_offsets[k] -= nx + 2;
	for (size_t i = 0; i < ny + 1; ++i)
    {
        for (size_t j = 0; j < nx + 1; ++j)
        {
            m_points[i * (nx + 1) + j] = Point(p1.x + hx * j, p1.y + hy * i);
        }
    }
    for (int i = 1; i < ny; ++i)
    {
        for (int k = (i + 1) * (nx + 1); k < size_p; ++k)
            m_offsets[k] -= 2;
    }
	for (size_t i = 0; i < ny; ++i)
	{
		for (size_t j = 0; j < nx; ++j)
		{
			Point temp[3];
			int nodes[3];
			nodes[0] = i * (nx + 1) + j;
			nodes[1] = i * (nx + 1) + j + 1;
			nodes[2] = (i + 1) * (nx + 1) + j + 1;
			temp[0] = m_points[nodes[0]];
			temp[1] = m_points[nodes[1]];
			temp[2] = m_points[nodes[2]];
            m_elems.push_back(new CFiniteElement < CTriangle, CTriangleBasis>{ nodes, temp, 1, 0 });

			nodes[0] = i * (nx + 1) + j;
			nodes[1] = (i + 1) * (nx + 1) + j;
			nodes[2] = (i + 1) * (nx + 1) + j + 1;
			temp[0] = m_points[nodes[0]];
			temp[1] = m_points[nodes[1]];
			temp[2] = m_points[nodes[2]];
            m_elems.push_back(new CFiniteElement < CTriangle, CTriangleBasis>{ nodes, temp, 1, 0 });
		}
	}
	Point temp[2];
	int nodes[2];
	for (size_t i = 0; i < nx; ++i)
	{
		temp[0] = m_points[i];
		temp[1] = m_points[i + 1];
		nodes[0] = i;
		nodes[1] = i + 1;
        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;
		m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 0 });
		m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * i);
		m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
		m_edges[m_edges.size() - 1]->ReverseNormal();
		nodes[0] = ny * (nx + 1) + i;
		nodes[1] = ny * (nx + 1) + i + 1;
		temp[0] = m_points[nodes[0]];
		temp[1] = m_points[nodes[1]];
        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;
		//m_edges[m_edges.size() - 1]->ReverseNormal();

		m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 1 });
		m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * (ny - 1) * nx + 2 * i + 1);
		m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
	}

	for (size_t i = 0; i < ny; ++i)
	{
		nodes[0] = i * (nx + 1);
		nodes[1] = (i + 1) * (nx + 1);
		temp[0] = m_points[nodes[0]];
		temp[1] = m_points[nodes[1]];
        m_bnds[nodes[0]] = 1;
        m_bnds[nodes[1]] = 1;
		m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, 2 });
		m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * i * nx + 1);
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
		m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * (nx - 1) + 2 * i * nx);
		m_edges[m_edges.size() - 1]->SetNeighbour(1, -1);
	}
	for (size_t i = 0; i < ny; ++i)
	{
		for (size_t j = 0; j < nx - 1; ++j)
		{
			nodes[0] = 1 + j + (nx + 1) * i;
			nodes[1] = 1 + j + (nx + 1) * (i + 1);
			temp[0] = m_points[nodes[0]];
			temp[1] = m_points[nodes[1]];
			m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
			m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * j + 2 * i * nx);
			m_edges[m_edges.size() - 1]->SetNeighbour(1, 2 * (j + 1) + 2 * i * nx + 1);
		}
	}

	for (size_t i = 0; i < ny - 1; ++i)
	{
		for (size_t j = 0; j < nx; ++j)
		{
			nodes[0] = nx + 1 + j + (nx + 1) * i;
			nodes[1] = nx + 2 + j + (nx + 1) * i;
			temp[0] = m_points[nodes[0]];
			temp[1] = m_points[nodes[1]];
			m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
			m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * j + 1 + 2 * i * nx);
			m_edges[m_edges.size() - 1]->SetNeighbour(1, 2 * j + 2 * (i + 1) * nx);
		}
	}
	
	for (size_t i = 0; i < ny; ++i)
	{
		for (size_t j = 0; j < nx; ++j)
		{
			nodes[0] = j + i * (nx + 1);
			nodes[1] = j + 1 + (i + 1) * (nx + 1);
			temp[0] = m_points[nodes[0]];
			temp[1] = m_points[nodes[1]];
			m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nodes, temp, 2, -1 });
			m_edges[m_edges.size() - 1]->SetNeighbour(0, 2 * j + 2 * i * nx);
			m_edges[m_edges.size() - 1]->SetNeighbour(1, 2 * j + 1 + 2 * i * nx);
		}
	}

    m_basepoints = m_points;
    const auto& sz_el = m_elems.size();
	const auto& sz_ed = m_edges.size();
    m_elemsbase.resize(sz_el);
	m_edgesbase.resize(sz_ed);
    int nods[3];
	Point points[3];
	for (int i = 0; i < sz_el; ++i)
	{
		nods[0] = m_elems[i]->GetNode(0);
		nods[1] = m_elems[i]->GetNode(1);
		nods[2] = m_elems[i]->GetNode(2);
		points[0] = m_points[nods[0]];
		points[1] = m_points[nods[1]];
		points[2] = m_points[nods[2]];
		m_elemsbase[i] = new CFiniteElement<CTriangle, CTriangleBasis>{ &nods[0], points, 1, m_elems[i]->GetType() };
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

const int CTriangularMesh::interpolate(const int node) const
{
    /*for (int i = 0; i < m_offsets.size(); ++i)
        cout << m_offsets[i] << endl;*/
    //if (node == 9)
        //cout << 9 << endl;
    if (m_bnds[node])
        return -1;
    return node + m_offsets[node];
}

const int CTriangularMesh::GetNumberOfINodes() const
{
    auto n = sqrt(m_points.size());
    return n * n - 4 * n + 4;
}

const int CTriangularMesh::refine_h()
{
	const int size_el = m_elems.size();
	int size_p = m_points.size();
	const int size_base = m_points.size();
	const int size_e = m_edges.size();
	vector<vector<int>>	elfs(size_el);
	vector<CElement<>*> n_elems;
	vector<CElement<>*> n_edges;
	vector<Point> n_points;
	vector<int> neis;
	vector<int> flags(size_el);
	int nodes[2];
	for (int i = 0; i < size_e; ++i)
	{
		const auto& edge = m_edges[i];
		vector<Point> points(2);
		points[0] = m_points[edge->GetNode(0)];
		points[1] = m_points[edge->GetNode(1)];
		Point temp;
		temp = Point{ points[0].x + (points[1].x - points[0].x) / 2, points[0].y + (points[1].y - points[0].y) / 2 };
		nodes[0] = edge->GetNode(0);
		nodes[1] = size_p;
		points[1] = temp;
		CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
		//edg1.SetNeighbour(0, edge->GetNeighbour(0));
		//edg1.SetNeighbour(1, edge->GetNeighbour(1));
		nodes[1] = edge->GetNode(1);
		nodes[0] = size_p;
		points[1] = m_points[edge->GetNode(1)];
		points[0] = temp;
		m_points.push_back(temp);
		CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
		//edg2.SetNeighbour(0, edge->GetNeighbour(0));
		//edg2.SetNeighbour(1, edge->GetNeighbour(1));
		++size_p;
		n_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg1 });
		n_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg2 });
		if(edge->GetNeighbour(0) > -1)
			elfs[edge->GetNeighbour(0)].push_back(i);
		if (edge->GetNeighbour(1) > -1)
			elfs[edge->GetNeighbour(1)].push_back(i);
	}
	const auto& ned = n_edges.size();

	for (int i = 0; i < size_el; ++i)
	{
		const auto& elem = m_elems[i];
		Point points[3];
		int nds[3];
		Point tp[3];
		int tnds[3];
		points[0] = m_points[elem->GetNode(0)];
		points[1] = m_points[elem->GetNode(1)];
		points[2] = m_points[elem->GetNode(2)];
		Point temp[3];
		temp[0] = Point{ points[0].x + (points[1].x - points[0].x) / 2, points[0].y + (points[1].y - points[0].y) / 2 };
		temp[1] = Point{ points[1].x + (points[2].x - points[1].x) / 2, points[1].y + (points[2].y - points[1].y) / 2 };
		temp[2] = Point{ points[0].x + (points[2].x - points[0].x) / 2, points[0].y + (points[2].y - points[0].y) / 2 };
		int j = size_base;
		for (; j < size_p; ++j)
			if (temp[0] == m_points[j])
				break;
		nds[0] = j;
		j = size_base;
		for (; j < size_p; ++j)
			if (temp[1] == m_points[j])
				break;
		nds[1] = j;
		j = size_base;
		for (; j < size_p; ++j)
			if (temp[2] == m_points[j])
				break;
		nds[2] = j;
		tp[0] = points[0];
		tp[1] = temp[0];
		tp[2] = temp[2];
		tnds[0] = elem->GetNode(0);
		tnds[1] = nds[0];
		tnds[2] = nds[2];
		n_elems.push_back(new CFiniteElement<CTriangle, CTriangleBasis>{ tnds, tp, 1, elem->GetType() });
		tp[0] = temp[0];
		tp[1] = points[1];
		tp[2] = temp[1];
		tnds[0] = nds[0];
		tnds[1] = elem->GetNode(1);
		tnds[2] = nds[1];
		n_elems.push_back(new CFiniteElement<CTriangle, CTriangleBasis>{ tnds, tp, 1, elem->GetType() });

		tp[0] = temp[1];
		tp[1] = points[2];
		tp[2] = temp[2];
		tnds[0] = nds[1];
		tnds[1] = elem->GetNode(2);
		tnds[2] = nds[2];
		n_elems.push_back(new CFiniteElement<CTriangle, CTriangleBasis>{ tnds, tp, 1, elem->GetType() });

		tp[0] = temp[0];
		tp[1] = temp[1];
		tp[2] = temp[2];
		tnds[0] = nds[0];
		tnds[1] = nds[1];
		tnds[2] = nds[2];
		n_elems.push_back(new CFiniteElement<CTriangle, CTriangleBasis>{ tnds, tp, 1, elem->GetType() });
	}
	auto scal = [&](const Point& p1, const Point& p2) {return p1.x*p2.x + p1.y*p2.y + p1.z * p2.z; };
	for (int i = 0; i < n_elems.size(); ++i)
	{
		int nods[3];
		int nds[3];
		Point points[3];
		int k = 0;
		const auto& elem = n_elems[i];
		nods[0] = elem->GetNode(0);
		nods[1] = elem->GetNode(1);
		nods[2] = elem->GetNode(2);
		nds[0] = nods[2];
		nds[1] = nods[0];
		points[0] = m_points[nds[0]];
		points[1] = m_points[nds[1]];
		points[2] = m_points[nods[1]];
		k = 0;
		CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nds, points, 2, -1} };
		for (; k < n_edges.size(); ++k)
		{
			if (edg1 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(n_edges[k]))
			{
				if (n_edges[k]->GetNeighbour(0) > -1)
					n_edges[k]->SetNeighbour(1, i);
				else
					n_edges[k]->SetNeighbour(0, i);
				break;
			}
		}
		if (k == n_edges.size())
		{
			n_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg1 });
			n_edges[n_edges.size() - 1]->SetNeighbour(0, i);
			if (scal(n_edges[n_edges.size() - 1]->GetNormal(), points[2]) > 0)
				n_edges[n_edges.size() - 1]->ReverseNormal();
		}
		nds[0] = nods[1];
		nds[1] = nods[2];
		points[0] = m_points[nds[0]];
		points[1] = m_points[nds[1]];
		points[2] = m_points[nods[0]];
		k = 0;
		CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nds, points, 2, -1} };
		for (; k < n_edges.size(); ++k)
		{
			if (edg2 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(n_edges[k]))
			{
				if (n_edges[k]->GetNeighbour(0) > -1)
					n_edges[k]->SetNeighbour(1, i);
				else
					n_edges[k]->SetNeighbour(0, i);
				break;
			}
		}
		if (k == n_edges.size())
		{
			n_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg2 });
			n_edges[n_edges.size() - 1]->SetNeighbour(0, i);
			if (scal(n_edges[n_edges.size() - 1]->GetNormal(), points[2]) > 0)
				n_edges[n_edges.size() - 1]->ReverseNormal();
		}

		points[0] = m_points[nods[0]];
		points[1] = m_points[nods[1]];
		points[2] = m_points[nods[2]];
		k = 0;
		CFiniteElement<CEdge, CEdgeLinearBasis, int> edg3{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nods, points, 2, -1} };

		for (; k < n_edges.size(); ++k)
		{
			if (edg3 == *static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(n_edges[k]))
			{
				if (n_edges[k]->GetNeighbour(0) > -1)
					n_edges[k]->SetNeighbour(1, i);
				else
					n_edges[k]->SetNeighbour(0, i);
				break;
			}
		}
		if (k == n_edges.size())
		{
			n_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ edg3 });
			n_edges[n_edges.size() - 1]->SetNeighbour(0, i);
			if (scal(n_edges[n_edges.size() - 1]->GetNormal(), points[2]) > 0)
				n_edges[n_edges.size() - 1]->ReverseNormal();
		}
	}

	auto sz(m_elems.size());
	for (auto i = 0; i < sz; ++i)
		delete m_elems[i];
	std::swap(m_elems, n_elems);
	vector<CElement<>*>().swap(n_elems);
	sz = m_edges.size();
	for (auto i = 0; i < sz; ++i)
		delete m_edges[i];
	std::swap(m_edges, n_edges);
	vector<CElement<>*>().swap(n_edges);
	vector<Point>().swap(m_basepoints);
	m_basepoints = m_points;
	if (m_elemsbase.size() > 0)
	{
		auto sz(m_elemsbase.size());
		for (auto i = 0; i < sz; ++i)
			delete m_elemsbase[i];
	}
	if (m_edgesbase.size() > 0)
	{
		auto sz(m_edgesbase.size());
		for (auto i = 0; i < sz; ++i)
			delete m_edgesbase[i];
	}
	const auto& sz_el = m_elems.size();
	const auto& sz_ed = m_edges.size();
	m_elemsbase.resize(sz_el);
	m_edgesbase.resize(sz_ed);
	int nods[3];
	Point points[3];
	for (int i = 0; i < sz_el; ++i)
	{
		//nods[0] = m_elems[i]->GetNode(0);
		//nods[1] = m_elems[i]->GetNode(1);
		//nods[2] = m_elems[i]->GetNode(2);
		//points[0] = m_points[nods[0]];
		//points[1] = m_points[nods[1]];
		//points[2] = m_points[nods[2]];
		//m_elemsbase[i] = new CFiniteElement<CTriangle, CTriangleBasis>{ nods, points, 3, m_elems[i]->GetType() };
		m_elemsbase[i] = m_elems[i]->Clone();
	}
	for (int i = 0; i < sz_ed; ++i)
	{
		//nods[0] = m_edges[i]->GetNode(0);
		//nods[1] = m_edges[i]->GetNode(1);
		//points[0] = m_points[nods[0]];
		//points[1] = m_points[nods[1]];
		//m_edgesbase[i] = new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ nods, points, 2, m_edges[i]->GetType() };
		m_edgesbase[i] = m_edges[i]->Clone();
	}

	/*for (int i = 0; i < ned; i += 2)
	{
		const auto& edge = n_edges[i];
		const int n1 = edge->GetNeighbour(0);
		const int n2 = edge->GetNeighbour(1);
		if (!flags[n1])
		{
			for (int j = i; j < ned; j += 2)
			{
				const auto& edge2 = n_edges[j];
				const int n11 = edge2->GetNeighbour(0);
				const int n22 = edge2->GetNeighbour(1);
				if (n11 == n1)
				{
					if (edge2->GetNode(0) == edge->GetNode(0))
					{
						Point points[3];
						int nds[3];
						points[0] = m_points[edge2->GetNode(0)];
					}
					flags[n1] = 1;
					break;
				}
			}
		}
	}*/
	return 0;
}

const int CTriangularMesh::refine_p()
{
	const auto order = m_order;//m_elems[0]->GetDoFs();
	const auto new_order = order + 1;
	m_order = new_order;

	/*if (m_elems.size() > 0)
	{
		auto sz(m_elems.size());
		for (auto i = 0; i < sz; ++i)
			delete m_elems[i];
	}
	if (m_edges.size() > 0)
	{
		auto sz(m_edges.size());
		for (auto i = 0; i < sz; ++i)
			delete m_edges[i];
	}

	const auto& sz_el = m_elems.size();
	const auto& sz_ed = m_edges.size();
	//m_elemsbase.resize(sz_el);
	//m_edgesbase.resize(sz_ed);
	int n_nods[3];
	Point n_points[3];
	for (int i = 0; i < sz_el; ++i)
	{
		n_nods[0] = m_elemsbase[i]->GetNode(0);
		n_nods[1] = m_elemsbase[i]->GetNode(1);
		n_nods[2] = m_elemsbase[i]->GetNode(2);
		n_points[0] = m_points[n_nods[0]];
		n_points[1] = m_points[n_nods[1]];
		n_points[2] = m_points[n_nods[2]];
		m_elems[i] = new CFiniteElement<CTriangle, CTriangleBasis>{ n_nods, n_points, 3, m_elemsbase[i]->GetType() };
	}
	for (int i = 0; i < sz_ed; ++i)
	{
		n_nods[0] = m_edgesbase[i]->GetNode(0);
		n_nods[1] = m_edgesbase[i]->GetNode(1);
		n_points[0] = m_points[n_nods[0]];
		n_points[1] = m_points[n_nods[1]];
		m_edges[i] = new CFiniteElement<CEdge, CEdgeLinearBasis, int>{ n_nods, n_points, 2, m_edgesbase[i]->GetType() };
	}*/

	//m_elems = m_elemsbase;
	//m_edges = m_edgesbase;
	const int size_el = m_elems.size();
	const int size_e = m_edges.size();
	vector<vector<int>>	elfs(size_el);
	vector<CElement<>*> n_elems;
	//vector<Point> n_points;
	vector<int> neis;
	vector<int> flags(size_el);
	int nodes[3];
	//m_points = m_basepoints;
	for (int i = 0; i < size_el; ++i)
		//for (int k = 0; k < order; ++k)
		m_elems[i]->IncreaseOrder();
	const auto s = 3 * (order)+(order - 1) * (order - 2) / 2;
	for (int i = 0; i < size_e; ++i)
	{
		const auto& edge = m_edges[i];
		vector<Point> points(2);
		points[0] = m_points[edge->GetNode(0)];
		points[1] = m_points[edge->GetNode(1)];
		vector<Point> temp(order);
		//for(int k = 0; k < order; ++k)
			m_edges[i]->IncreaseOrder();
		//for (int j = 0; j < order; ++j)
		{
			const auto szm = m_points.size();
			temp[0] = Point{ points[0].x + (points[1].x - points[0].x) / new_order, points[0].y + (points[1].y - points[0].y) / new_order };
			m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
			const auto& nk = m_edges[i]->GetNeighbour(0);
			if (nk > -1)
			{
				nodes[0] = m_elems[nk]->GetNode(0);
				nodes[1] = m_elems[nk]->GetNode(1);
				//nodes[2] = m_elems[nk]->GetNode(2);
				CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
				nodes[0] = m_elems[nk]->GetNode(1);
				nodes[1] = m_elems[nk]->GetNode(2);
				CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
				nodes[0] = m_elems[nk]->GetNode(2);
				nodes[1] = m_elems[nk]->GetNode(0);
				CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
				int number = 0;
				if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
					number = 0;
				else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
					number = 1;
				else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
					number = 2;
				switch (number)
				{
				case 0:
					m_elems[nk]->SetNode(s, szm);
					break;
				case 1:
					m_elems[nk]->SetNode(s + 1, szm);
					break;
				case 2:
					m_elems[nk]->SetNode(s + 2, szm);
					break;
				default:
					break;
				}
			}
			const auto& ne = m_edges[i]->GetNeighbour(1);
			if(ne > - 1)
			{
				nodes[0] = m_elems[ne]->GetNode(0);
				nodes[1] = m_elems[ne]->GetNode(1);
				//nodes[2] = m_elems[nk]->GetNode(2);
				CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
				nodes[0] = m_elems[ne]->GetNode(1);
				nodes[1] = m_elems[ne]->GetNode(2);
				CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
				nodes[0] = m_elems[ne]->GetNode(2);
				nodes[1] = m_elems[ne]->GetNode(0);
				CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
				int number = 0;
				if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
					number = 0;
				else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
					number = 1;
				else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
					number = 2;
				switch (number)
				{
				case 0:
					m_elems[ne]->SetNode(s, szm);
					break;
				case 1:
					m_elems[ne]->SetNode(s + 1, szm);
					break;
				case 2:
					m_elems[ne]->SetNode(s + 2, szm);
					break;
				default:
					break;
				}
			}
		}
	}
	const int size_base = m_basepoints.size();
	const int size_p = m_points.size();
	const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
	//const int centerpr = 0.5 + double(new_order - 3) / 2.;
	//const int center = 0.5 + double(new_order - 2) / 2.;
	const int center = s2 - s - 3;
    cout << center << endl;
    if (center == 1)
    {
        for (int i = 0; i < size_el; ++i)
        {
            Point pts[3];
            for (int k = 0; k < center; ++k)
            {
                m_elems[i]->SetNode(s + 3 + k, m_points.size());
                pts[0] = m_points[m_elems[i]->GetNode(0)];
                pts[2] = m_points[m_elems[i]->GetNode(1)];
                pts[1] = m_points[m_elems[i]->GetNode(2)];

                Point temp;
                temp.x = (pts[0].x + pts[1].x + pts[2].x) / 3;
                temp.y = (pts[0].y + pts[1].y + pts[2].y) / 3;
                m_points.push_back(temp);
            }
        }
    }
    else
    {
        for (int i = 0; i < size_el; ++i)
        {
            /*const auto& elem = m_elems[i];
            Point points[3];
            int nds[3];
            Point tp[3];
            int tnds[3];
            points[0] = m_points[elem->GetNode(0)];
            points[1] = m_points[elem->GetNode(1)];
            points[2] = m_points[elem->GetNode(2)];
            vector<vector<Point>> temp(3);
            const auto dofs = m_elems[i]->GetDoFs();
            //for (int k = 0; k < order; ++k)
            {
                temp[0].push_back( Point{ points[0].x + (points[1].x - points[0].x) / new_order, points[0].y + (points[1].y - points[0].y) / new_order });
                temp[1].push_back(Point{ points[1].x + (points[2].x - points[1].x) / new_order, points[1].y + (points[2].y - points[1].y) / new_order });
                temp[2].push_back(Point{ points[0].x + (points[2].x - points[0].x) / new_order, points[0].y + (points[2].y - points[0].y) / new_order });
                m_elems[i]->IncreaseOrder();
            }*/
            double r1, r2;
            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_real_distribution<double> dist(0., 1.);
            Point pts[3];
            for (int k = 0; k < center; ++k)
            {
                m_elems[i]->SetNode(s + 3 + k, m_points.size());
                r1 = dist(mt);
                r2 = dist(mt);
                pts[0] = m_points[m_elems[i]->GetNode(0)];
                pts[2] = m_points[m_elems[i]->GetNode(1)];
                pts[1] = m_points[m_elems[i]->GetNode(2)];

                Point temp;
                temp.x = (1 - GaussTriangle::m_tra[k] - GaussTriangle::m_trb[k]) * pts[0].x + GaussTriangle::m_tra[k] * pts[1].x + GaussTriangle::m_trb[k] * pts[2].x;
                temp.y = (1 - GaussTriangle::m_tra[k] - GaussTriangle::m_trb[k]) * pts[0].y + GaussTriangle::m_tra[k] * pts[1].y + GaussTriangle::m_trb[k] * pts[2].y;

                temp = Point((1 - sqrt(r1))*pts[0].x + sqrt(r1)*(1 - r2)*pts[1].x + sqrt(r1)*r2*pts[2].x, (1 - sqrt(r1))*pts[0].y + sqrt(r1)*(1 - r2)*pts[1].y + sqrt(r1)*r2*pts[2].y);
                m_points.push_back(temp);
            }
            /*const auto& dofs2 = m_elems[i]->GetDoFs();
            for(int k = 0; k < 3; ++k)
                //for (int m = dofs; m < dofs2; ++m)
                //for(int m = 0; m < 3; ++m)
                    for (int j = size_base; j < size_p; ++j)
                        if (temp[k][0] == m_points[j])
                        {
                            m_elems[i]->SetNode(dofs + k, j);
                            break;
                        }*/
        }
    }
	return 0;
}


const int CTriangularMesh::set4thOrder()
{
    set4thNodes_1();
    set4thNodes_2();

    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    const int size_el = m_elems.size();
    const int size_e = m_edges.size();
    vector<vector<int>>	elfs(size_el);
    vector<CElement<>*> n_elems;
    n_elems.resize(size_el);
    vector<int> neis;
    vector<int> flags(size_el);
    int nodes[3];
    for (int i = 0; i < size_el; ++i)
        m_elems[i]->IncreaseOrder();
    const auto s = 3 * (order)+(order - 1) * (order - 2) / 2 - 1;
    for (int i = 0; i < size_e; ++i)
    {
        const auto& edge = m_edges[i];
        vector<Point> points(2);
        points[0] = m_points[edge->GetNode(0)];
        points[1] = m_points[edge->GetNode(1)];
        vector<Point> temp(order);
        //for(int k = 0; k < order; ++k)
            m_edges[i]->IncreaseOrder();
        //for (int j = 0; j < order; ++j)
        {
            const auto szm = m_points.size();
            temp[0] = Point{ points[0].x + 3. * (points[1].x - points[0].x) / 4, points[0].y + 3. *(points[1].y - points[0].y) / 4 };
            m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
            const auto& nk = m_edges[i]->GetNeighbour(0);
            if (nk > -1)
            {
                nodes[0] = m_elems[nk]->GetNode(0);
                nodes[1] = m_elems[nk]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(1);
                nodes[1] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(2);
                nodes[1] = m_elems[nk]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[nk]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[nk]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[nk]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
            const auto& ne = m_edges[i]->GetNeighbour(1);
            if(ne > - 1)
            {
                nodes[0] = m_elems[ne]->GetNode(0);
                nodes[1] = m_elems[ne]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(1);
                nodes[1] = m_elems[ne]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(2);
                nodes[1] = m_elems[ne]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[ne]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[ne]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[ne]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
        }
    }
    const int size_base = m_basepoints.size();
    const int size_p = m_points.size();
    const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
    const int center = s2 - s - 3;
    Point pts[3];
    Point temp;
    Point mid1, mid2, mid3;
    for (int i = 0; i < size_el; ++i)
    {
        pts[0] = m_points[m_elems[i]->GetNode(0)];
        pts[1] = m_points[m_elems[i]->GetNode(1)];
        pts[2] = m_points[m_elems[i]->GetNode(2)];
        m_elems[i]->SetNode(s + 3, m_points.size());
        mid1.x = (pts[0].x + pts[1].x) / 2;
        mid1.y = (pts[0].y + pts[1].y) / 2;

        mid2.x = (pts[0].x + pts[2].x) / 2;
        mid2.y = (pts[0].y + pts[2].y) / 2;

        temp.x = (mid1.x + mid2.x) / 2;
        temp.y = (mid1.y + mid2.y) / 2;
        m_points.push_back(temp);

        m_elems[i]->SetNode(s + 4, m_points.size());
        mid3.x = (pts[1].x + pts[2].x) / 2;
        mid3.y = (pts[1].y + pts[2].y) / 2;

        temp.x = (mid1.x + mid3.x) / 2;
        temp.y = (mid1.y + mid3.y) / 2;
        m_points.push_back(temp);
        m_elems[i]->SetNode(s + 5, m_points.size());
        temp.x = (mid3.x + mid2.x) / 2;
        temp.y = (mid3.y + mid2.y) / 2;
        m_points.push_back(temp);
    }
    int nods[20];
    Point points[20];
    for (int i = 0; i < size_el; ++i)
    {
        nods[0] = m_elems[i]->GetNode(0);
        nods[1] = m_elems[i]->GetNode(1);
        nods[2] = m_elems[i]->GetNode(2);
        nods[3] = m_elems[i]->GetNode(3);
        nods[4] = m_elems[i]->GetNode(4);
        nods[5] = m_elems[i]->GetNode(5);
        nods[6] = m_elems[i]->GetNode(6);
        nods[7] = m_elems[i]->GetNode(7);
        nods[8] = m_elems[i]->GetNode(8);
        nods[9] = m_elems[i]->GetNode(9);
        nods[10] = m_elems[i]->GetNode(10);
        nods[11] = m_elems[i]->GetNode(11);
        nods[12] = m_elems[i]->GetNode(12);
        nods[13] = m_elems[i]->GetNode(13);
        nods[14] = m_elems[i]->GetNode(14);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        points[3] = m_points[nods[3]];
        points[4] = m_points[nods[4]];
        points[5] = m_points[nods[5]];
        points[6] = m_points[nods[6]];
        points[7] = m_points[nods[7]];
        points[8] = m_points[nods[8]];
        points[9] = m_points[nods[9]];
        points[10] = m_points[nods[10]];
        points[11] = m_points[nods[11]];
        points[12] = m_points[nods[12]];
        points[13] = m_points[nods[13]];
        points[14] = m_points[nods[14]];
        n_elems[i] = new CFiniteElement<CTriangle, CTriangleLagrangeBasis>{ &nods[0], points, 4, m_elems[i]->GetType() };
        //n_elems[i]->IncreaseOrder();
        //n_elems[i]->IncreaseOrder();
        //n_elems[i]->IncreaseOrder();
    }

    if(m_elems.size() > 0)
    {
        auto sz(m_elems.size());
        for(auto i = 0; i < sz; ++i)
            delete m_elems[i];
    }
    for(int i = 0; i < size_el; ++i)
    {
        nods[0] = n_elems[i]->GetNode(0);
        nods[1] = n_elems[i]->GetNode(1);
        nods[2] = n_elems[i]->GetNode(2);
        nods[3] = n_elems[i]->GetNode(3);
        nods[4] = n_elems[i]->GetNode(4);
        nods[5] = n_elems[i]->GetNode(5);
        nods[6] = n_elems[i]->GetNode(6);
        nods[7] = n_elems[i]->GetNode(7);
        nods[8] = n_elems[i]->GetNode(8);
        nods[9] = n_elems[i]->GetNode(9);
        nods[10] = n_elems[i]->GetNode(10);
        nods[11] = n_elems[i]->GetNode(11);
        nods[12] = n_elems[i]->GetNode(12);
        nods[13] = n_elems[i]->GetNode(13);
        nods[14] = n_elems[i]->GetNode(14);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        points[3] = m_points[nods[3]];
        points[4] = m_points[nods[4]];
        points[5] = m_points[nods[5]];
        points[6] = m_points[nods[6]];
        points[7] = m_points[nods[7]];
        points[8] = m_points[nods[8]];
        points[9] = m_points[nods[9]];
        points[10] = m_points[nods[10]];
        points[11] = m_points[nods[11]];
        points[12] = m_points[nods[12]];
        points[13] = m_points[nods[13]];
        points[14] = m_points[nods[14]];
        m_elems[i] = new CFiniteElement<CTriangle, CTriangleLagrangeBasis>{ &nods[0], points, 4, n_elems[i]->GetType() };
        //m_elems[i]->IncreaseOrder();
        //m_elems[i]->IncreaseOrder();
        //m_elems[i]->IncreaseOrder();
    }
    if (n_elems.size() > 0)
    {
        auto sz(n_elems.size());
        for (auto i = 0; i < sz; ++i)
            delete n_elems[i];
        vector<CElement<>*>().swap(n_elems);
    }
    return 0;
}


const int CTriangularMesh::set2ndOrder()
{
    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    const int size_el = m_elems.size();
    const int size_e = m_edges.size();
    vector<vector<int>>	elfs(size_el);
    vector<CElement<>*> n_elems;
    n_elems.resize(size_el);
    vector<int> neis;
    vector<int> flags(size_el);
    int nodes[3];
    for (int i = 0; i < size_el; ++i)
        m_elems[i]->IncreaseOrder();
    const auto s = 3 * (order)+(order - 1) * (order - 2) / 2;
    for (int i = 0; i < size_e; ++i)
    {
        const auto& edge = m_edges[i];
        vector<Point> points(2);
        points[0] = m_points[edge->GetNode(0)];
        points[1] = m_points[edge->GetNode(1)];
        vector<Point> temp(order);
        //for(int k = 0; k < order; ++k)
            m_edges[i]->IncreaseOrder();
        //for (int j = 0; j < order; ++j)
        {
            const auto szm = m_points.size();
            temp[0] = Point{ points[0].x + (points[1].x - points[0].x) / new_order, points[0].y + (points[1].y - points[0].y) / new_order };
            m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
            const auto& nk = m_edges[i]->GetNeighbour(0);
            if (nk > -1)
            {
                nodes[0] = m_elems[nk]->GetNode(0);
                nodes[1] = m_elems[nk]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(1);
                nodes[1] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(2);
                nodes[1] = m_elems[nk]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[nk]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[nk]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[nk]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
            const auto& ne = m_edges[i]->GetNeighbour(1);
            if(ne > - 1)
            {
                nodes[0] = m_elems[ne]->GetNode(0);
                nodes[1] = m_elems[ne]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(1);
                nodes[1] = m_elems[ne]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(2);
                nodes[1] = m_elems[ne]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[ne]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[ne]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[ne]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
        }
    }
    const int size_base = m_basepoints.size();
    const int size_p = m_points.size();
    const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
    const int center = s2 - s - 3;
    for (int i = 0; i < size_el; ++i)
    {
        for (int k = 0; k < center; ++k)
        {
            m_elems[i]->SetNode(s + 3 + k, m_points.size());
            m_points.push_back(Point(0, 0, 0));
        }
    }
    int nods[20];
    Point points[20];
    for (int i = 0; i < size_el; ++i)
    {
        nods[0] = m_elems[i]->GetNode(0);
        nods[1] = m_elems[i]->GetNode(1);
        nods[2] = m_elems[i]->GetNode(2);
        nods[3] = m_elems[i]->GetNode(3);
        nods[4] = m_elems[i]->GetNode(4);
        nods[5] = m_elems[i]->GetNode(5);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        points[3] = m_points[nods[3]];
        points[4] = m_points[nods[4]];
        points[5] = m_points[nods[5]];
        n_elems[i] = new CFiniteElement<CTriangle, CTriangleLagrangeBasis>{ &nods[0], points, 2, m_elems[i]->GetType() };
        //n_elems[i]->IncreaseOrder();
        //n_elems[i]->IncreaseOrder();
        //n_elems[i]->IncreaseOrder();
    }

    if(m_elems.size() > 0)
    {
        auto sz(m_elems.size());
        for(auto i = 0; i < sz; ++i)
            delete m_elems[i];
    }
    for(int i = 0; i < size_el; ++i)
    {
        nods[0] = n_elems[i]->GetNode(0);
        nods[1] = n_elems[i]->GetNode(1);
        nods[2] = n_elems[i]->GetNode(2);
        nods[3] = n_elems[i]->GetNode(3);
        nods[4] = n_elems[i]->GetNode(4);
        nods[5] = n_elems[i]->GetNode(5);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        points[3] = m_points[nods[3]];
        points[4] = m_points[nods[4]];
        points[5] = m_points[nods[5]];
        m_elems[i] = new CFiniteElement<CTriangle, CTriangleLagrangeBasis>{ &nods[0], points, 2, n_elems[i]->GetType() };
        //m_elems[i]->IncreaseOrder();
        //m_elems[i]->IncreaseOrder();
        //m_elems[i]->IncreaseOrder();
    }
    if (n_elems.size() > 0)
    {
        auto sz(n_elems.size());
        for (auto i = 0; i < sz; ++i)
            delete n_elems[i];
        vector<CElement<>*>().swap(n_elems);
    }
    return 0;
}
const unsigned int CTriangularMesh::GetNumberOfNodes() const
{
	return (unsigned int)m_points.size();
}

const int CTriangularMesh::set3rdOrder()
{
    //refine_p();
    set3rdNodes();
    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    const int size_el = m_elems.size();
    const int size_e = m_edges.size();
    vector<vector<int>>	elfs(size_el);
    vector<CElement<>*> n_elems;
    n_elems.resize(size_el);
    vector<int> neis;
    vector<int> flags(size_el);
    int nodes[3];
    for (int i = 0; i < size_el; ++i)
        m_elems[i]->IncreaseOrder();
    const auto s = 3 * (order)+(order - 1) * (order - 2) / 2;
    for (int i = 0; i < size_e; ++i)
    {
        const auto& edge = m_edges[i];
        vector<Point> points(2);
        points[0] = m_points[edge->GetNode(0)];
        points[1] = m_points[edge->GetNode(1)];
        vector<Point> temp(order);
        //for(int k = 0; k < order; ++k)
            m_edges[i]->IncreaseOrder();
        //for (int j = 0; j < order; ++j)
        {
            const auto szm = m_points.size();
            temp[0] = Point{ points[0].x + 2. * (points[1].x - points[0].x) / new_order, points[0].y + 2. * (points[1].y - points[0].y) / new_order };
            m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
            const auto& nk = m_edges[i]->GetNeighbour(0);
            if (nk > -1)
            {
                nodes[0] = m_elems[nk]->GetNode(0);
                nodes[1] = m_elems[nk]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(1);
                nodes[1] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(2);
                nodes[1] = m_elems[nk]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[nk]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[nk]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[nk]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
            const auto& ne = m_edges[i]->GetNeighbour(1);
            if(ne > - 1)
            {
                nodes[0] = m_elems[ne]->GetNode(0);
                nodes[1] = m_elems[ne]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(1);
                nodes[1] = m_elems[ne]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(2);
                nodes[1] = m_elems[ne]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[ne]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[ne]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[ne]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
        }
    }
    const int size_base = m_basepoints.size();
    const int size_p = m_points.size();
    const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
    const int center = s2 - s - 3;
    for (int i = 0; i < size_el; ++i)
    {
        for (int k = 0; k < center; ++k)
        {
            m_elems[i]->SetNode(s + 3 + k, m_points.size());
            //m_points.push_back({0, 0, 0});
            m_points.push_back(Point((m_points[m_elems[i]->GetNode(0)].x + m_points[m_elems[i]->GetNode(1)].x + m_points[m_elems[i]->GetNode(2)].x) / 3, (m_points[m_elems[i]->GetNode(0)].y + m_points[m_elems[i]->GetNode(1)].y + m_points[m_elems[i]->GetNode(2)].y) / 3, 0));
        }
    }
    int nods[20];
    Point points[20];
    cout << "here" << endl;
    for (int i = 0; i < size_el; ++i)
    {
        nods[0] = m_elems[i]->GetNode(0);
        nods[1] = m_elems[i]->GetNode(1);
        nods[2] = m_elems[i]->GetNode(2);
        nods[3] = m_elems[i]->GetNode(3);
        nods[4] = m_elems[i]->GetNode(4);
        nods[5] = m_elems[i]->GetNode(5);
        nods[6] = m_elems[i]->GetNode(6);
        nods[7] = m_elems[i]->GetNode(7);
        nods[8] = m_elems[i]->GetNode(8);
        nods[9] = m_elems[i]->GetNode(9);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        points[3] = m_points[nods[3]];
        points[4] = m_points[nods[4]];
        points[5] = m_points[nods[5]];
        points[6] = m_points[nods[6]];
        points[7] = m_points[nods[7]];
        points[8] = m_points[nods[8]];
        points[9] = m_points[nods[9]];
        n_elems[i] = new CFiniteElement<CTriangle, CTriangleLagrangeBasis>{ &nods[0], points, 3, m_elems[i]->GetType() };
        //n_elems[i]->IncreaseOrder();
        //n_elems[i]->IncreaseOrder();
        //n_elems[i]->IncreaseOrder();
    }

    if(m_elems.size() > 0)
    {
        auto sz(m_elems.size());
        for(auto i = 0; i < sz; ++i)
            delete m_elems[i];
    }
    cout << "here2" << endl;
    for(int i = 0; i < size_el; ++i)
    {
        nods[0] = n_elems[i]->GetNode(0);
        nods[1] = n_elems[i]->GetNode(1);
        nods[2] = n_elems[i]->GetNode(2);
        nods[3] = n_elems[i]->GetNode(3);
        nods[4] = n_elems[i]->GetNode(4);
        nods[5] = n_elems[i]->GetNode(5);
        nods[6] = n_elems[i]->GetNode(6);
        nods[7] = n_elems[i]->GetNode(7);
        nods[8] = n_elems[i]->GetNode(8);
        nods[9] = n_elems[i]->GetNode(9);
        points[0] = m_points[nods[0]];
        points[1] = m_points[nods[1]];
        points[2] = m_points[nods[2]];
        points[3] = m_points[nods[3]];
        points[4] = m_points[nods[4]];
        points[5] = m_points[nods[5]];
        points[6] = m_points[nods[6]];
        points[7] = m_points[nods[7]];
        points[8] = m_points[nods[8]];
        points[9] = m_points[nods[9]];
        m_elems[i] = new CFiniteElement<CTriangle, CTriangleLagrangeBasis>{ &nods[0], points, 3, n_elems[i]->GetType() };
        //m_elems[i]->IncreaseOrder();
        //m_elems[i]->IncreaseOrder();
        //m_elems[i]->IncreaseOrder();
    }
    cout << "done" << endl;
    if (n_elems.size() > 0)
    {
        auto sz(n_elems.size());
        for (auto i = 0; i < sz; ++i)
            delete n_elems[i];
        vector<CElement<>*>().swap(n_elems);
    }
    return 0;
}

const Mesh::Point CTriangularMesh::GetNode(const unsigned int n) const
{
	return m_points[n];
}

void CTriangularMesh::set3rdNodes()
{
    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    const int size_el = m_elems.size();
    const int size_e = m_edges.size();
    vector<vector<int>>	elfs(size_el);
    vector<CElement<>*> n_elems;
    //vector<Point> n_points;
    vector<int> neis;
    vector<int> flags(size_el);
    int nodes[3];
    for (int i = 0; i < size_el; ++i)
        m_elems[i]->IncreaseOrder();
    const auto s = 3 * (order)+(order - 1) * (order - 2) / 2;
    for (int i = 0; i < size_e; ++i)
    {
        const auto& edge = m_edges[i];
        vector<Point> points(2);
        points[0] = m_points[edge->GetNode(0)];
        points[1] = m_points[edge->GetNode(1)];
        vector<Point> temp(order);
            m_edges[i]->IncreaseOrder();
        //for (int j = 0; j < order; ++j)
        {
            const auto szm = m_points.size();
            temp[0] = Point{ points[0].x + (points[1].x - points[0].x) / 3, points[0].y + (points[1].y - points[0].y) / 3 };
            m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
            const auto& nk = m_edges[i]->GetNeighbour(0);
            if (nk > -1)
            {
                nodes[0] = m_elems[nk]->GetNode(0);
                nodes[1] = m_elems[nk]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(1);
                nodes[1] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(2);
                nodes[1] = m_elems[nk]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[nk]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[nk]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[nk]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
            const auto& ne = m_edges[i]->GetNeighbour(1);
            if(ne > - 1)
            {
                nodes[0] = m_elems[ne]->GetNode(0);
                nodes[1] = m_elems[ne]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(1);
                nodes[1] = m_elems[ne]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(2);
                nodes[1] = m_elems[ne]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[ne]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[ne]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[ne]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
        }
    }
    const int size_base = m_basepoints.size();
    const int size_p = m_points.size();
    const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
    const int center = s2 - s - 3;
    for (int i = 0; i < size_el; ++i)
    {
        for (int k = 0; k < center; ++k)
        {
            m_elems[i]->SetNode(s + 3 + k, m_points.size());
            m_points.push_back(Point(0, 0, 0));
        }
    }
}


const unsigned int CTriangularMesh::GetNumberOfElements() const
{
	return (unsigned int)m_elems.size();
}

void CTriangularMesh::set4thNodes_1()
{
    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    const int size_el = m_elems.size();
    const int size_e = m_edges.size();
    vector<vector<int>>	elfs(size_el);
    vector<CElement<>*> n_elems;
    //vector<Point> n_points;
    vector<int> neis;
    vector<int> flags(size_el);
    int nodes[3];
    for (int i = 0; i < size_el; ++i)
        m_elems[i]->IncreaseOrder();
    const auto s = 3 * (order)+(order - 1) * (order - 2) / 2;
    for (int i = 0; i < size_e; ++i)
    {
        const auto& edge = m_edges[i];
        vector<Point> points(2);
        points[0] = m_points[edge->GetNode(0)];
        points[1] = m_points[edge->GetNode(1)];
        vector<Point> temp(order);
            m_edges[i]->IncreaseOrder();
        //for (int j = 0; j < order; ++j)
        {
            const auto szm = m_points.size();
            temp[0] = Point{ points[0].x + (points[1].x - points[0].x) / 4, points[0].y + (points[1].y - points[0].y) / 4 };
            m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
            const auto& nk = m_edges[i]->GetNeighbour(0);
            if (nk > -1)
            {
                nodes[0] = m_elems[nk]->GetNode(0);
                nodes[1] = m_elems[nk]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(1);
                nodes[1] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(2);
                nodes[1] = m_elems[nk]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[nk]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[nk]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[nk]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
            const auto& ne = m_edges[i]->GetNeighbour(1);
            if(ne > - 1)
            {
                nodes[0] = m_elems[ne]->GetNode(0);
                nodes[1] = m_elems[ne]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(1);
                nodes[1] = m_elems[ne]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(2);
                nodes[1] = m_elems[ne]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[ne]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[ne]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[ne]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
        }
    }
    const int size_base = m_basepoints.size();
    const int size_p = m_points.size();
    const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
    const int center = s2 - s - 3;
    /*for (int i = 0; i < size_el; ++i)
    {
        for (int k = 0; k < center; ++k)
        {
            m_elems[i]->SetNode(s + 3 + k, m_points.size());
            m_points.push_back(Point(0, 0, 0));
        }
    }*/
}

void CTriangularMesh::set4thNodes_2()
{
    const auto order = m_order;//m_elems[0]->GetDoFs();
    const auto new_order = order + 1;
    m_order = new_order;
    const int size_el = m_elems.size();
    const int size_e = m_edges.size();
    vector<vector<int>>	elfs(size_el);
    vector<CElement<>*> n_elems;
    //vector<Point> n_points;
    vector<int> neis;
    vector<int> flags(size_el);
    int nodes[3];
    for (int i = 0; i < size_el; ++i)
        m_elems[i]->IncreaseOrder();
    const auto s = 3 * (order)+(order - 1) * (order - 2) / 2;
    for (int i = 0; i < size_e; ++i)
    {
        const auto& edge = m_edges[i];
        vector<Point> points(2);
        points[0] = m_points[edge->GetNode(0)];
        points[1] = m_points[edge->GetNode(1)];
        vector<Point> temp(order);
            m_edges[i]->IncreaseOrder();
        //for (int j = 0; j < order; ++j)
        {
            const auto szm = m_points.size();
            temp[0] = Point{ points[0].x + (points[1].x - points[0].x) / 2, points[0].y + (points[1].y - points[0].y) / 2 };
            m_edges[i]->SetNode(1 + order, szm);
            m_points.push_back(temp[0]);
            const auto& nk = m_edges[i]->GetNeighbour(0);
            if (nk > -1)
            {
                nodes[0] = m_elems[nk]->GetNode(0);
                nodes[1] = m_elems[nk]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(1);
                nodes[1] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[nk]->GetNode(2);
                nodes[1] = m_elems[nk]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[nk]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[nk]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[nk]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
            const auto& ne = m_edges[i]->GetNeighbour(1);
            if(ne > - 1)
            {
                nodes[0] = m_elems[ne]->GetNode(0);
                nodes[1] = m_elems[ne]->GetNode(1);
                //nodes[2] = m_elems[nk]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg0{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(1);
                nodes[1] = m_elems[ne]->GetNode(2);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                nodes[0] = m_elems[ne]->GetNode(2);
                nodes[1] = m_elems[ne]->GetNode(0);
                CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{ CFiniteElement<CEdge, CEdgeLinearBasis, int>{nodes, &points[0], 2, edge->GetType()} };
                int number = 0;
                if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg0)
                    number = 0;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg1)
                    number = 1;
                else if (*static_cast<CFiniteElement<CEdge, CEdgeLinearBasis, int>*>(m_edges[i]) == edg2)
                    number = 2;
                switch (number)
                {
                case 0:
                    m_elems[ne]->SetNode(s, szm);
                    break;
                case 1:
                    m_elems[ne]->SetNode(s + 1, szm);
                    break;
                case 2:
                    m_elems[ne]->SetNode(s + 2, szm);
                    break;
                default:
                    break;
                }
            }
        }
    }
    const int size_base = m_basepoints.size();
    const int size_p = m_points.size();
    const auto s2 = 3 * (order + 1)+(order) * (order - 1) / 2;
    const int center = s2 - s - 3;
}
const int CTriangularMesh::FindElement(const Point& test) const
{
	Point p[3];
	double volElem, mvol;
	for (unsigned int i = 0; i < m_elems.size(); ++i)
	{
		mvol = 0;
		//CFiniteElement<TriangleQuad, QuadBasisTr> elem(m_fems[i]);
		const auto& elem = m_elems[i];
		//CFiniteElement<TriangleQuad, LinearBasisTr> elem(m_fems[i]);
		p[0] = m_points[elem->GetNode(0)];
		p[1] = m_points[elem->GetNode(1)];
		p[2] = m_points[elem->GetNode(2)];
		volElem = elem->GetMeasure();
		mvol = fabs(CompSquare(p[0], p[1], test));
		mvol += fabs(CompSquare(p[1], p[2], test));
		mvol += fabs(CompSquare(p[0], p[2], test));
		if (fabs(mvol - volElem) < 1e-13)
		{
			return i;
		}
	}
	return -1;
}

const CElement<>* CTriangularMesh::GetElement(const unsigned int n) const
{
	if(n < m_elems.size())
		return m_elems[n];
	return nullptr;
}

const CElement<>* CTriangularMesh::GetBoundary(const unsigned int n) const
{
	if(n < m_edges.size())
		return m_edges[n];
	return nullptr;
}

const double CTriangularMesh::getSolution(const unsigned int element, const unsigned int node) const
{
	return 0.0;
}

const int CTriangularMesh::updateSolution(const unsigned int element, const unsigned int node, const double value)
{
	return 0;
}

const std::vector<double> CTriangularMesh::getSolution() const
{
	return std::vector<double>();
}

const int CTriangularMesh::updateSolution(const std::vector<double>&)
{
	return 0;
}

const int CTriangularMesh::updateSolution(const unsigned int element, const unsigned int node, CSolution * value)
{
	return 0;
}

const double CTriangularMesh::getParameter(Parameters param, const unsigned int l, const Point & p) const
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

const double CTriangularMesh::getParameter(Parameters param, const unsigned int l, const int i) const
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

const int CTriangularMesh::setParameter(Parameters param, const double, const unsigned int)
{
	return 0;
}

const int CTriangularMesh::setParameter(const CParameter& p, const unsigned int type)
{
	m_params[type] = p;
	return 0;
}

const int CTriangularMesh::updateSolution(const unsigned int node, const double value)
{
	return 0;
}


const unsigned int CTriangularMesh::GetNumberOfBoundaries() const
{
	return (unsigned int)m_edges.size();
}
