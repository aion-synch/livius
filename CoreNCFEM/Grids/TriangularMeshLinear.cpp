#include <stdio.h>
#include "TriangularMeshLinear.h"
#include "../FiniteElements/TriangleLinear.h"
#include "../FiniteElements/Edge.h"
#include <iostream>
#include <algorithm>
using namespace std;
using namespace corenc;
using namespace Mesh;

CTriangularMeshLinear::CTriangularMeshLinear()
{
	
}

CTriangularMeshLinear::~CTriangularMeshLinear()
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
}

CTriangularMeshLinear::CTriangularMeshLinear(const string& file_name)
{
    cout << file_name << endl;
    cout << "Reading the grid from the file..." << endl;
    ifstream in{file_name};
    if (!in.is_open())
    {
        cout << "Error! Couldn't load the grid\n";
        return;
    }
	CFiniteElement<CTriangleLinear, CTriangleLinearBasis, int> ele{};
	CFiniteElement<CEdge, CEdgeLinearBasis, int> edge{};
	CEdge edg;
	CEdgeLinearBasis edglb;
    float temp;
    int num, tp;
    int nums;
    int nods[3], nds[3];
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
            cout << "\nError: size is less than 1" << endl;
            throw;
        }
        m_points.resize(num);
        for (int i = 0; i < num; ++i)
            in >> tmp >> m_points[i].x >> m_points[i].y >> tmp;
        in >> tmp >> tmp;
        in >> nums;
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
					points[0] = m_points[nods[0]];
					points[1] = m_points[nods[1]];
					m_edges.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, int>(nods, points, 2, tp));
                    break;
				case 2:
				{
					in >> tp >> tp >> tp >> nods[0] >> nods[1] >> nods[2];
					--nods[0]; --nods[1]; --nods[2];
					nds[0] = nods[2];
					nds[1] = nods[0];
					points[0] = m_points[nds[0]];
					points[1] = m_points[nds[1]];
					points[2] = m_points[nods[1]];
					k = 0;
					CFiniteElement<CEdge, CEdgeLinearBasis, int> edg1{CFiniteElement<CEdge, CEdgeLinearBasis, int>{nds, points, 2, -1}};
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
					CFiniteElement<CEdge, CEdgeLinearBasis, int> edg2{CFiniteElement<CEdge, CEdgeLinearBasis, int>{nds, points, 2, -1}};
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
					CFiniteElement<CEdge, CEdgeLinearBasis, int> edg3{CFiniteElement<CEdge, CEdgeLinearBasis, int>{nods, points, 2, -1}};
					
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
					m_elems.push_back(new CFiniteElement<CTriangleLinear, CTriangleLinearBasis, int>{nods, points, 3, tp});
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
        cout << "Sorry, error: ";
        cout << e.what() << endl;
        exit(0);
    }
}

corenc::Mesh::CTriangularMeshLinear::CTriangularMeshLinear(const CTriangularMeshLinear &)
{
}

const int CTriangularMeshLinear::refine_h()
{
	const int size_el = m_elems.size();
	int size_p = m_points.size();
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
		int j = 0;
		for (; j < size_p; ++j)
			if (temp[0] == m_points[j])
				break;
		nds[0] = j;
		j = 0;
		for (; j < size_p; ++j)
			if (temp[1] == m_points[j])
				break;
		nds[1] = j;
		j = 0;
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
		n_elems.push_back(new CFiniteElement<CTriangleLinear, CTriangleLinearBasis, int>{ tnds, tp, 3, elem->GetType() });

		tp[0] = temp[0];
		tp[1] = points[1];
		tp[2] = temp[1];
		tnds[0] = nds[0];
		tnds[1] = elem->GetNode(1);
		tnds[2] = nds[1];
		n_elems.push_back(new CFiniteElement<CTriangleLinear, CTriangleLinearBasis, int>{ tnds, tp, 3, elem->GetType() });

		tp[0] = temp[1];
		tp[1] = points[2];
		tp[2] = temp[2];
		tnds[0] = nds[1];
		tnds[1] = elem->GetNode(2);
		tnds[2] = nds[2];
		n_elems.push_back(new CFiniteElement<CTriangleLinear, CTriangleLinearBasis, int>{ tnds, tp, 3, elem->GetType() });

		tp[0] = temp[0];
		tp[1] = temp[1];
		tp[2] = temp[2];
		tnds[0] = nds[0];
		tnds[1] = nds[1];
		tnds[2] = nds[2];
		n_elems.push_back(new CFiniteElement<CTriangleLinear, CTriangleLinearBasis, int>{ tnds, tp, 3, elem->GetType() });
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
	m_elems.swap(n_elems);
	vector<CElement<>*>().swap(n_elems);
	sz = m_edges.size();
	for (auto i = 0; i < sz; ++i)
		delete m_edges[i];
	m_edges.swap(n_edges);
	vector<CElement<>*>().swap(n_edges);
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

const unsigned int CTriangularMeshLinear::GetNumberOfNodes() const
{
	return (unsigned int)m_points.size();
}

const Mesh::Point CTriangularMeshLinear::GetNode(const unsigned int n) const
{
	return m_points[n];
}
const unsigned int CTriangularMeshLinear::GetNumberOfElements() const
{
	return (unsigned int)m_elems.size();
}

const int CTriangularMeshLinear::FindElement(const Mesh::Point&) const
{
	return -1;
}

const CElement<>* CTriangularMeshLinear::GetElement(const unsigned int n) const
{
	if(n < m_elems.size())
		return m_elems[n];
	return nullptr;
}

const CElement<>* CTriangularMeshLinear::GetBoundary(const unsigned int n) const
{
	if(n < m_edges.size())
		return m_edges[n];
	return nullptr;
}

const double CTriangularMeshLinear::getSolution(const unsigned int element, const unsigned int node) const
{
	return 0.0;
}

const int CTriangularMeshLinear::updateSolution(const unsigned int element, const unsigned int node, const double value)
{
	return 0;
}

const std::vector<double> CTriangularMeshLinear::getSolution() const
{
	return std::vector<double>();
}

const int CTriangularMeshLinear::updateSolution(const std::vector<double>&)
{
	return 0;
}

const int CTriangularMeshLinear::updateSolution(const unsigned int element, const unsigned int node, CSolution * value)
{
	return 0;
}

const double CTriangularMeshLinear::getParameter(Parameters param, const unsigned int l, const Point & p) const
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

const double CTriangularMeshLinear::getParameter(Parameters param, const unsigned int l, const int i) const
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

const int CTriangularMeshLinear::setParameter(Parameters param, const double, const unsigned int)
{
	return 0;
}

const int CTriangularMeshLinear::setParameter(const CParameter& p, const unsigned int type)
{
	m_params[type] = p;
	return 0;
}

const int CTriangularMeshLinear::updateSolution(const unsigned int node, const double value)
{
	return 0;
}


const unsigned int CTriangularMeshLinear::GetNumberOfBoundaries() const
{
	return (unsigned int)m_edges.size();
}
