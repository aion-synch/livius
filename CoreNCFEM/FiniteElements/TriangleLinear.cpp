#include "TriangleLinear.h"
#include <iostream>
using namespace corenc;
using namespace Mesh;
using namespace std;

CTriangleLinear::CTriangleLinear()
{
    m_nodes[0] = m_nodes[1] = m_nodes[2] = m_edges[0] = m_edges[1] = m_edges[2] = -1;
}

CTriangleLinear::CTriangleLinear(const int n1, const int n2, const int n3)
{
    m_nodes[0] = n1;
    m_nodes[1] = n2;
    m_nodes[2] = n3;
    m_edges[0] = m_edges[1] = m_edges[2] = -1;
}

CTriangleLinear::CTriangleLinear(const int n1, const int n2, const int n3, const int e1, const int e2, const int e3)
{
    m_nodes[0] = n1;
    m_nodes[1] = n2;
    m_nodes[2] = n3;
    m_edges[0] = e1;
    m_edges[1] = e2;
    m_edges[2] = e3;
}

CTriangleLinear::CTriangleLinear(const int* nodes)
{
    m_nodes[0] = nodes[0];
    m_nodes[1] = nodes[1];
    m_nodes[2] = nodes[2];
}

CTriangleLinear::CTriangleLinear(const int* nodes, const int* edges)
{
    m_nodes[0] = nodes[0];
    m_nodes[1] = nodes[1];
    m_nodes[2] = nodes[2];
    m_edges[0] = edges[0];
    m_edges[1] = edges[1];
    m_edges[2] = edges[2];
}

CTriangleLinear::CTriangleLinear(const CTriangleLinear& t)
{
    m_nodes[0] = t.m_nodes[0];
    m_nodes[1] = t.m_nodes[1];
    m_nodes[2] = t.m_nodes[2];
    m_edges[0] = t.m_edges[0];
    m_edges[1] = t.m_edges[1];
    m_edges[2] = t.m_edges[2];
}

const int CTriangleLinear::GetNode(const int n) const
{
    if(n < 3)
        return m_nodes[n];
    cout << "GetNode(" << n << "): Wrong number of the node." << endl;
    return -1;
}

const int CTriangleLinear::GetNode(const NODES& node) const
{
	if (node == NODES::FIRST)
		return m_nodes[0];
	return m_nodes[1];
}

const int CTriangleLinear::GetEdge(const int n) const
{
    if(n < 3)
        return m_edges[n];
    cout << "GetEdge(" << n << "): Wrong number of the edge." << endl;
    return -1;
}

const int CTriangleLinear::GetFacet(const int) const
{
    return -1;
}

const int CTriangleLinear::GetNumberOfNodes() const
{
    return 3;
}

const int CTriangleLinear::GetNumberOfEdges() const
{
    return 3;
}

const int CTriangleLinear::GetNumberOfFacets() const
{
    return 1;
}

void CTriangleLinear::SetNode(const int n, const int node)
{
    if(n < 3)
        m_nodes[n] = node;
    else
        cout << "SetNode(" << n << ", " << node << "): Wrong number of the node." << endl;
}

void CTriangleLinear::SetEdge(const int k, const int edge)
{
    if(k < 3)
        m_edges[k] = edge;
    else
        cout << "SetEdge(" << k << ", " << edge << "): Wrong number of the edge." << endl;
}

void CTriangleLinear::SetFacet(const int, const int)
{
    cout << "You tried." << endl;
}

const double CTriangleLinear::Integrate(const std::function<const double (const Point &)> &f, const std::vector<Point>& v) const
{
	double py, pz, det;
	det = (v[1].y - v[0].y)*(v[2].z - v[0].z) - (v[2].y - v[0].y)*(v[1].z - v[0].z);
	det *= det;
	py = (v[2].x - v[0].x)*(v[1].z - v[0].z) - (v[1].x - v[0].x)*(v[2].z - v[0].z);
	py *= py;
	pz = (v[1].x - v[0].x)*(v[2].y - v[0].y) - (v[2].x - v[0].x)*(v[1].y - v[0].y);
	pz *= pz;
	det = sqrt(det + py + pz);
	Point temp;
	double integral = 0;
	for (int i = 0; i < 7; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian();
	}
	return det * integral;
}

const Point CTriangleLinear::Integrate(const std::function<const Point (const Point &)> &f, const std::vector<Point>& v) const
{
	Point temp;
	Point integral{0,0,0};
	for (int i = 0; i < 7; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian();
	}
	return integral;
}

const vector<double> CTriangleLinear::Integrate(const std::function<const std::vector<double> (const Point &)> &f, const std::vector<Point> &v) const
{
	Point temp;
	auto size {f({0,0,0}).size()};
	vector<double> integral;
	integral.resize(size);
	for (int i = 0; i < 7; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		auto fx {f(temp)};
		for(int j = 0; j < size; ++j)
			integral[j] = GaussTriangle::m_trw[i] * fx[j]*temp.Jacobian();
	}
	return integral;
}

CTriangleLinearBasis::CTriangleLinearBasis()
{
    m_det = 0;
    m_alpha[0][0] = m_alpha[0][1] = m_alpha[0][2] =
    m_alpha[1][0] = m_alpha[1][1] = m_alpha[1][2] =
    m_alpha[2][0] = m_alpha[2][1] = m_alpha[2][2] = 0;
	m_normal = Point(0, 0, 0);
}

CTriangleLinearBasis::CTriangleLinearBasis(const Point& p1, const Point& p2, const Point& p3)
{
    compD(p1, p2, p3);
	compAlpha(p1, p2, p3);
	compNormal(p1, p2, p3);
}

CTriangleLinearBasis::CTriangleLinearBasis(const Point* p)
{
	compD(p[0], p[1], p[2]);
	compAlpha(p[0], p[1], p[2]);
	compNormal(p[0], p[1], p[2]);
}

const int CTriangleLinearBasis::GetNumberOfShapeFunctions() const
{
	return m_number;
}

const Point CTriangleLinearBasis::GetNormal() const
{
	return m_normal;
}

void CTriangleLinearBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}

const double CTriangleLinearBasis::GetShapeFunction(const int k, const Point &p) const
{
	//const DForm<0> lp = [&](const Point& p1){return (m_alpha[k][0] + m_alpha[k][1] * p.x + m_alpha[k][2] * p.y);};
	if(k < 3)
		return (m_alpha[k][0] + m_alpha[k][1] * p.x + m_alpha[k][2] * p.y);
	cout << "CTriangleLinearBasis::GetShapeFunction(" << k << "): is out of range." << endl;
	return -1;
}

const Point CTriangleLinearBasis::GetGradShapeFunction(const int k, const Point &) const
{
	if(k < 3)
		return Point{m_alpha[k][1], m_alpha[k][2], 0};
	cout << "CTriangleLinearBasis::GetGradShapeFunction(" << k << "): is out of range." << endl;
	return Point(-1,-1,-1);
}
//const function<const double(const Point&)> CTriangleLinearBasis::GetShapeFunction(const int k) const
//{
	//return [&](const Point& p){return m_alpha[k][0] + m_alpha[k][1] * p.x + m_alpha[k][2] * p.y;};
//}

CTriangleLinearBasis::CTriangleLinearBasis(const CTriangleLinearBasis& t)
{
    m_normal = t.m_normal;
    m_det = t.m_det;
	m_alpha[0][0] = t.m_alpha[0][0];
	m_alpha[0][1] = t.m_alpha[0][1];
	m_alpha[0][2] = t.m_alpha[0][2];
	
	m_alpha[1][0] = t.m_alpha[1][0];
	m_alpha[1][1] = t.m_alpha[1][1];
	m_alpha[1][2] = t.m_alpha[1][2];
	
	m_alpha[2][0] = t.m_alpha[2][0];
	m_alpha[2][1] = t.m_alpha[2][1];
	m_alpha[2][2] = t.m_alpha[2][2];
}

void CTriangleLinearBasis::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
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

void CTriangleLinearBasis::compAlpha(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
{
	m_alpha[0][0] = (p2.x * p3.y - p2.y * p3.x) / m_det;
	m_alpha[1][0] = -(p1.x * p3.y - p1.y * p3.x) / m_det;
	m_alpha[2][0] = (p1.x * p2.y - p1.y * p2.x) / m_det;
	
	m_alpha[0][1] = -(p3.y - p2.y) / m_det;
	m_alpha[1][1] = (p3.y - p1.y) / m_det;
	m_alpha[2][1] = -(p2.y - p1.y) / m_det;
	
	m_alpha[0][2] = (p3.x - p2.x) / m_det;
	m_alpha[1][2] = -(p3.x - p1.x) / m_det;
	m_alpha[2][2] = (p2.x - p1.x) / m_det;
}

void CTriangleLinearBasis::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
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

const double CTriangleLinearBasis::GetValue(const Point& p) const
{
	return 0.;
}
//const int CTriangleLinearBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const CFESolution CTriangleLinearBasis::GetValue(const int n) const
//{
//	return m_w[n];
//}



































