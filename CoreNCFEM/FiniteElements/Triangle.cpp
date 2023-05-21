#include "Triangle.h"
#include <iostream>
#include <algorithm>
#include <random>
#include "../../CoreNCA/Matrix.h"
#include "../../CoreNCA/MatrixSkyline.h"
using namespace corenc;
using namespace Mesh;
using namespace std;
using namespace Algebra;

CTriangle::CTriangle()
{
    m_edges[0] = m_edges[1] = m_edges[2] = -1;
	m_nodes.resize(3);
	m_order = 0;
	m_number = 0;
}

CTriangle::CTriangle(const int n1, const int n2, const int n3, const int order)
{
	m_nodes.resize(3);
    m_nodes[0] = n1;
    m_nodes[1] = n2;
    m_nodes[2] = n3;
	m_number = order;
    m_edges[0] = m_edges[1] = m_edges[2] = -1;
	SetOrder();
}

CTriangle::CTriangle(const int n1, const int n2, const int n3, const int e1, const int e2, const int e3, const int order)
{
	m_nodes.resize(3);
    m_nodes[0] = n1;
    m_nodes[1] = n2;
    m_nodes[2] = n3;
    m_edges[0] = e1;
    m_edges[1] = e2;
    m_edges[2] = e3;
	m_number = order;
	SetOrder();
}

CTriangle::CTriangle(const int* nodes, const int order)
{
	m_number = order;
	SetOrder();
	m_nodes.resize(m_order);
	for (int i = 0; i < m_order; ++i)
		m_nodes[i] = nodes[i];
}

CTriangle::CTriangle(const int* nodes, const int* edges, const int order)
{
	m_number = order;
	SetOrder();
	m_nodes.resize(m_order);
	for (int i = 0; i < m_order; ++i)
		m_nodes[i] = nodes[i];
    m_edges[0] = edges[0];
    m_edges[1] = edges[1];
    m_edges[2] = edges[2];
	m_number = 1;
}

CTriangle::CTriangle(const CTriangle& t)
{
    m_nodes = t.m_nodes;
	m_order = t.m_order;
    m_edges[0] = t.m_edges[0];
    m_edges[1] = t.m_edges[1];
    m_edges[2] = t.m_edges[2];
	m_number = t.m_number;
}

const int CTriangle::GetNode(const int n) const
{
    if(n < m_nodes.size())
        return m_nodes[n];
    cout << "GetNode(" << n << "): Wrong number of the node." << endl;
    return -1;
}

const int CTriangle::GetNode(const NODES& node) const
{
	if (node == NODES::FIRST)
		return m_nodes[0];
	return m_nodes[1];
}

const int CTriangle::GetEdge(const int n) const
{
    if(n < 3)
        return m_edges[n];
    cout << "GetEdge(" << n << "): Wrong number of the edge." << endl;
    return -1;
}

const int CTriangle::GetFacet(const int) const
{
    return -1;
}

const int CTriangle::GetNumberOfNodes() const
{
    return (int)m_nodes.size();
}

const int CTriangle::GetNumberOfEdges() const
{
    return 3;
}

const int CTriangle::GetNumberOfFacets() const
{
    return 1;
}

void CTriangle::SetNode(const int n, const int node)
{
    if(n < m_order)
        m_nodes[n] = node;
    else
        cout << "SetNode(" << n << ", " << node << "): Wrong number of the node." << endl;
}

const int CTriangle::IncreaseOrder()
{
	++m_number;
	//const auto pr = m_order;
	//m_order = 3 * (m_number) + (m_number - 1) * (m_number - 2) / 2;
	SetOrder();
	m_nodes.resize(m_order);
	//for (int i = pr; i < m_order; ++i)
		//m_nodes[i] = -1;
	return 0;
}

void CTriangle::SetOrder()
{
	m_order = 3 * (m_number)+(m_number - 1) * (m_number - 2) / 2;
}

void CTriangle::SetEdge(const int k, const int edge)
{
    if(k < 3)
        m_edges[k] = edge;
    else
        cout << "SetEdge(" << k << ", " << edge << "): Wrong number of the edge." << endl;
}

void CTriangle::SetFacet(const int, const int)
{
    cout << "You tried." << endl;
}

const double CTriangle::Integrate(const std::function<const double (const Point &)> &f, const std::vector<Point>& v) const
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
	for (int i = 0; i < GaussTriangle::m_order; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian();
	}
	// don't look 
	/*Point tp[3];
	tp[0].x = v[0].x + (v[1].x - v[0].x) / 2;
	tp[0].y = v[0].y + (v[1].y - v[0].y) / 2;

	tp[1].x = v[1].x + (v[2].x - v[1].x) / 2;
	tp[1].y = v[1].y + (v[2].y - v[1].y) / 2;

	tp[2].x = v[2].x + (v[0].x - v[2].x) / 2;
	tp[2].y = v[2].y + (v[0].y - v[2].y) / 2;

	for (int i = 0; i < GaussTriangle::m_order; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * tp[0].x + GaussTriangle::m_trb[i] * tp[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * tp[0].y + GaussTriangle::m_trb[i] * tp[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * tp[0].z + GaussTriangle::m_trb[i] * tp[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian() / 4;
	}

	for (int i = 0; i < GaussTriangle::m_order; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * tp[1].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * tp[1].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * tp[1].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian() / 4;
	}

	for (int i = 0; i < GaussTriangle::m_order; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[2].x + GaussTriangle::m_tra[i] * tp[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[2].y + GaussTriangle::m_tra[i] * tp[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[2].z + GaussTriangle::m_tra[i] * tp[1].z + GaussTriangle::m_trb[i] * v[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian() / 4;
	}

	for (int i = 0; i < GaussTriangle::m_order; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[0].x + GaussTriangle::m_tra[i] * tp[1].x + GaussTriangle::m_trb[i] * tp[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[0].y + GaussTriangle::m_tra[i] * tp[1].y + GaussTriangle::m_trb[i] * tp[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * tp[0].z + GaussTriangle::m_tra[i] * tp[1].z + GaussTriangle::m_trb[i] * tp[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian() / 4;
	}*/

	return 0.5 * det * integral;
}

const Point CTriangle::Integrate(const std::function<const Point (const Point &)> &f, const std::vector<Point>& v) const
{
	Point temp;
	Point integral{0,0,0};
	for (int i = 0; i < 16; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian();
	}
	return integral;
}

const vector<double> CTriangle::Integrate(const std::function<const std::vector<double> (const Point &)> &f, const std::vector<Point> &v) const
{
	Point temp;
	auto size {f({0,0,0}).size()};
	vector<double> integral;
	integral.resize(size);
	for (int i = 0; i < 16; ++i)
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

CTriangleBasis::CTriangleBasis()
{
    m_det = 0;
    m_alpha[0][0] = m_alpha[0][1] = m_alpha[0][2] =
    m_alpha[1][0] = m_alpha[1][1] = m_alpha[1][2] =
    m_alpha[2][0] = m_alpha[2][1] = m_alpha[2][2] = 0;
	m_sp = 0;
	m_s = 0;
	m_1dorder = 0;
	m_order = 0;
	createS();
	m_normal = Point(0, 0, 0);
	m_number = 0;
	//m_all.resize(1);
}
CTriangleBasis::CTriangleBasis(const Point& p1, const Point& p2, const Point& p3, const int order)
{
    compD(p1, p2, p3);
	compAlpha(p1, p2, p3);
	compNormal(p1, p2, p3);
	m_order = order;
	createS();
	m_1dorder = 1;
	m_number = 3;
	//m_all.resize(1);
}

CTriangleBasis::CTriangleBasis(const Point* p, const int order)
{
	compD(p[0], p[1], p[2]);
	compAlpha(p[0], p[1], p[2]);
	compNormal(p[0], p[1], p[2]);
	m_order = order;
	createS();
	m_1dorder = 1;


}

const int CTriangleBasis::createS()
{
	if (m_order < 2)
	{
		m_s = 3;
		m_sp = 0;
		return 0;
	}
	m_sp = std::max(int(3 * (m_order - 1) + (m_order - 2) * (m_order - 3) / 2), 0);
	m_s = 3 * (m_order) + (m_order - 1) * (m_order - 2) / 2;
	return 0;
}

const int CTriangleBasis::GetNumberOfShapeFunctions() const
{
	return m_s;
}

const Point CTriangleBasis::GetNormal() const
{
	return m_normal;
}

void CTriangleBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}
const double CTriangleBasis::m_L(const int k, const Point &p) const
{
	return (m_alpha[k][0] + m_alpha[k][1] * p.x + m_alpha[k][2] * p.y);
}

const double CTriangleBasis::m_xi(const int k, const Point &p) const
{
	return std::pow(p.x, (k - 2));
}

const double CTriangleBasis::GetShapeFunction(const int k, const Point &p) const
{
	if (k < 3)
		return m_L(k, p);
	switch (k)
	{
		// quadratic
        case 3: return m_L(0, p) * m_L(1, p);
        case 4: return m_L(1, p) * m_L(2, p);
        case 5: return m_L(0, p) * m_L(2, p);
            // cubic
        case 6: return m_L(0, p) * m_L(1, p) * (m_L(0, p) - m_L(1, p));
        case 7: return m_L(1, p) * m_L(2, p) * (m_L(1, p) - m_L(2, p));
        case 8: return m_L(0, p) * m_L(2, p) * (m_L(0, p) - m_L(2, p));
        case 9: return m_L(0, p) * m_L(1, p) * m_L(2, p);
        case 10: return m_L(0, p) * m_L(1, p) * (m_L(0, p) - m_L(1, p)) * (m_L(0, p) - m_L(1, p));
        case 11: return m_L(1, p) * m_L(2, p) * (m_L(1, p) - m_L(2, p)) * (m_L(1, p) - m_L(2, p));
        case 12: return m_L(0, p) * m_L(2, p) * (m_L(0, p) - m_L(2, p)) * (m_L(0, p) - m_L(2, p));
        case 13: return m_L(0, p) * m_L(1, p) * m_L(2, p) * (2 * m_L(2, p) - 1);
        case 14: return m_L(0, p) * m_L(1, p) * m_L(2, p) * (m_L(0, p) - m_L(1, p));
        default:
            //return 0.;
            break;
	}
	const auto s = 3 * (m_1dorder - 1) + (m_1dorder - 2) * (m_1dorder - 3) / 2;
	// ki = s + k
	int i = 0;
	for (i = 0; i < m_1dorder - 1; ++i)
	{
		if (k < m_all[i])
			break;
	}
	int ki = k - m_all[i - 1];
	const auto d = m_1dorder - 2;
	switch (ki)
	{
	case 0: return m_L(0, p) * m_L(1, p) * std::pow(m_L(0, p) - m_L(1, p), d);
	case 1: return m_L(1, p) * m_L(2, p) * std::pow(m_L(1, p) - m_L(2, p), d);
	case 2: return m_L(0, p) * m_L(2, p) * std::pow(m_L(0, p) - m_L(2, p), d);
	//case 3: return m_L(0, p) * m_L(1, p) * m_L(2, p) * std::pow(2 * m_L(2, p) - 1, m_1dorder - ki);
	default: break;
	}
	return m_L(0, p) * m_L(1, p) * m_L(2, p) * std::pow(m_L(0, p) - m_L(1, p), ki - 3) * std::pow(2 * m_L(2, p) - 1, m_1dorder - ki);
}

const Point CTriangleBasis::GetGradShapeFunction(const int k, const Point & p) const
{
	if(k < 3)
		return Point{m_alpha[k][1], m_alpha[k][2], 0};
    double temp = 0;
	switch (k)
	{
		// quadratic
	case 3: return Point(m_alpha[0][1] * m_L(1, p) + m_alpha[1][1] * m_L(0, p), m_alpha[0][2] * m_L(1, p) + m_alpha[1][2] * m_L(0, p), 0);
	case 4: return Point(m_alpha[1][1] * m_L(2, p) + m_alpha[2][1] * m_L(1, p), m_alpha[1][2] * m_L(2, p) + m_alpha[2][2] * m_L(1, p), 0);
	case 5: return Point(m_alpha[0][1] * m_L(2, p) + m_alpha[2][1] * m_L(0, p), m_alpha[0][2] * m_L(2, p) + m_alpha[2][2] * m_L(0, p), 0);
		// cubic
	case 6: return Point(	m_L(1, p) * (m_L(0, p) - m_L(1, p)) * m_alpha[0][1] + m_L(0, p) * (m_L(0, p) - m_L(1, p)) * m_alpha[1][1] + m_L(0, p) * m_L(1, p) * (m_alpha[0][1] - m_alpha[1][1]),
							m_L(1, p) * (m_L(0, p) - m_L(1, p)) * m_alpha[0][2] + m_L(0, p) * (m_L(0, p) - m_L(1, p)) * m_alpha[1][2] + m_L(0, p) * m_L(1, p) * (m_alpha[0][2] - m_alpha[1][2]), 0);

	case 7: return Point(	m_L(2, p) * (m_L(1, p) - m_L(2, p)) * m_alpha[1][1] + m_L(1, p) * (m_L(1, p) - m_L(2, p)) * m_alpha[2][1] + m_L(1, p) * m_L(2, p) * (m_alpha[1][1] - m_alpha[2][1]),
							m_L(2, p) * (m_L(1, p) - m_L(2, p)) * m_alpha[1][2] + m_L(1, p) * (m_L(1, p) - m_L(2, p)) * m_alpha[2][2] + m_L(1, p) * m_L(2, p) * (m_alpha[1][2] - m_alpha[2][2]), 0);

    case 8: return Point(	m_L(2, p) * (m_L(0, p) - m_L(2, p)) * m_alpha[0][1] + m_L(0, p) * (m_L(0, p) - m_L(2, p)) * m_alpha[2][1] + m_L(0, p) * m_L(2, p) * (m_alpha[0][1] - m_alpha[2][1]),
                            m_L(2, p) * (m_L(0, p) - m_L(2, p)) * m_alpha[0][2] + m_L(0, p) * (m_L(0, p) - m_L(2, p)) * m_alpha[2][2] + m_L(0, p) * m_L(2, p) * (m_alpha[0][2] - m_alpha[2][2]), 0);

	case 9: return Point(	m_L(1, p) * m_L(2, p) * m_alpha[0][1] + m_L(0, p) * m_L(2, p) * m_alpha[1][1] + m_L(0, p) * m_L(1, p) * m_alpha[2][1],
							m_L(1, p) * m_L(2, p) * m_alpha[0][2] + m_L(0, p) * m_L(2, p) * m_alpha[1][2] + m_L(0, p) * m_L(1, p) * m_alpha[2][2], 0);

    case 10: temp = m_L(0, p) - m_L(1, p);
        return Point(m_alpha[0][1] * m_L(1, p) * temp * temp + m_alpha[1][1] * m_L(0, p) * temp * temp + 2 * (m_alpha[0][1] - m_alpha[1][1]) * temp * m_L(0, p) * m_L(1, p),
                m_alpha[0][2] * m_L(1, p) * temp * temp + m_alpha[1][2] * m_L(0, p) * temp * temp + 2 * (m_alpha[0][2] - m_alpha[1][2]) * temp * m_L(0, p) * m_L(1, p));
    case 11: temp = m_L(1, p) - m_L(2, p);
        return Point(m_alpha[1][1] * m_L(2, p) * temp * temp + m_alpha[2][1] * m_L(1, p) * temp * temp + 2 * (m_alpha[1][1] - m_alpha[2][1]) * temp * m_L(1, p) * m_L(2, p),
                m_alpha[1][2] * m_L(2, p) * temp * temp + m_alpha[2][2] * m_L(1, p) * temp * temp + 2 * (m_alpha[1][2] - m_alpha[2][2]) * temp * m_L(1, p) * m_L(2, p));
    case 12: temp = m_L(0, p) - m_L(2, p);
        return Point(m_alpha[2][1] * m_L(0, p) * temp * temp + m_alpha[0][1] * m_L(2, p) * temp * temp + 2 * (m_alpha[0][1] - m_alpha[2][1]) * temp * m_L(0, p) * m_L(2, p),
                m_alpha[2][2] * m_L(0, p) * temp * temp + m_alpha[0][2] * m_L(2, p) * temp * temp + 2 * (m_alpha[0][2] - m_alpha[2][2]) * temp * m_L(0, p) * m_L(2, p));
    case 13: temp = 2 * m_L(2, p) - 1;
        return Point(m_alpha[0][1] * m_L(1, p) * m_L(2, p) * temp + m_alpha[1][1] * m_L(0, p) * m_L(2, p) * temp + m_alpha[2][1] * m_L(0, p) * m_L(1, p) * temp + 2 * m_alpha[2][1] * m_L(0, p) * m_L(1, p) * m_L(2, p),
                m_alpha[0][2] * m_L(1, p) * m_L(2, p) * temp + m_alpha[1][2] * m_L(0, p) * m_L(2, p) * temp + m_alpha[2][2] * m_L(0, p) * m_L(1, p) * temp + 2 * m_alpha[2][2] * m_L(0, p) * m_L(1, p) * m_L(2, p));
    case 14: temp = m_L(0, p) - m_L(1, p);
        return Point(m_alpha[0][1] * m_L(1, p) * m_L(2, p) * temp + m_alpha[1][1] * m_L(0, p) * m_L(2, p) * temp + m_alpha[2][1] * m_L(0, p) * m_L(1, p) * temp + (m_alpha[0][1] - m_alpha[1][1]) * m_L(0, p) * m_L(1, p) * m_L(2, p),
                m_alpha[0][2] * m_L(1, p) * m_L(2, p) * temp + m_alpha[1][2] * m_L(0, p) * m_L(2, p) * temp + m_alpha[2][2] * m_L(0, p) * m_L(1, p) * temp + (m_alpha[0][2] - m_alpha[1][2]) * m_L(0, p) * m_L(1, p) * m_L(2, p));
	default: 
		break;
	}
	const auto s = 3 * (m_1dorder - 1) + (m_1dorder - 2) * (m_1dorder - 3) / 2;
	const auto sp = 3 * (m_1dorder) + (m_1dorder - 1) * (m_1dorder - 2) / 2;
	// ki = s + k
	int i = 0;
	for (i = 0; i < m_1dorder - 1; ++i)
	{
		if (k < m_all[i])
			break;
	}
	int ki = k - m_all[i - 1];
	double val, valp;
    const auto d = m_order - 2;
	switch (ki)
	{
	case 0: 
		valp = m_L(0, p) - m_L(1, p);
		val = std::pow(valp, d - 1);
		valp = val * valp;
		return Point(	m_L(1, p) * valp * m_alpha[0][1] + m_L(0, p) * valp * m_alpha[1][1] + d * m_L(0, p) * m_L(1, p) * (m_alpha[0][1] - m_alpha[1][1]) * val,
						m_L(1, p) * valp * m_alpha[0][2] + m_L(0, p) * valp * m_alpha[1][2] + d * m_L(0, p) * m_L(1, p) * (m_alpha[0][2] - m_alpha[1][2]) * val);
	case 1:
		valp = m_L(1, p) - m_L(2, p);
		val = std::pow(valp, d - 1);
		valp = val * valp;
		return Point(	m_L(2, p) * valp * m_alpha[1][1] + m_L(1, p) * valp * m_alpha[2][1] + d * m_L(1, p) * m_L(2, p) * (m_alpha[1][1] - m_alpha[2][1]) * val,
						m_L(2, p) * valp * m_alpha[1][2] + m_L(1, p) * valp * m_alpha[2][2] + d * m_L(1, p) * m_L(2, p) * (m_alpha[1][2] - m_alpha[2][2]) * val);
	case 2:
		valp = m_L(0, p) - m_L(2, p);
		val = std::pow(valp, d - 1);
		valp = val * valp;
		return Point(	m_L(2, p) * valp * m_alpha[0][1] + m_L(0, p) * valp * m_alpha[2][1] + d * m_L(0, p) * m_L(2, p) * (m_alpha[0][1] - m_alpha[2][1]) * val,
						m_L(2, p) * valp * m_alpha[0][2] + m_L(0, p) * valp * m_alpha[2][2] + d * m_L(0, p) * m_L(2, p) * (m_alpha[0][2] - m_alpha[2][2]) * val);
	default: break;
	}
	const int m = ki - 3;
	const int n = m_1dorder - ki;
	val = std::pow(m_L(0, p) - m_L(1, p), m);
	valp = std::pow(2 * m_L(2, p) - 1, n);
	double val2;// = std::pow(m_L(0, p) - m_L(1, p), m - 1);
	double valp2;// = std::pow(2 * m_L(2, p) - 1, n - 1);
	if (m == 0)
	{
		valp2 = std::pow(2 * m_L(2, p) - 1, n - 1);
		return Point(	valp * (m_L(1, p) * m_L(2, p) * m_alpha[0][1] + m_L(0, p) * m_L(2, p) * m_alpha[1][1] + m_L(0, p) * m_L(1, p) * m_alpha[2][1]) +
						m_L(0, p) * m_L(1, p) * m_L(2, p) * (2 * n * valp2 * m_alpha[2][1]),
						valp * (m_L(1, p) * m_L(2, p) * m_alpha[0][2] + m_L(0, p) * m_L(2, p) * m_alpha[1][2] + m_L(0, p) * m_L(1, p) * m_alpha[2][2]) +
						m_L(0, p) * m_L(1, p) * m_L(2, p) * (2 * n * valp2 * m_alpha[2][2]));
	}
	if (n == 0)
	{
		val2 = std::pow(m_L(0, p) - m_L(1, p), m - 1);
		return Point(	val * (m_L(1, p) * m_L(2, p) * m_alpha[0][1] + m_L(0, p) * m_L(2, p) * m_alpha[1][1] + m_L(0, p) * m_L(1, p) * m_alpha[2][1]) +
						m_L(0, p) * m_L(1, p) * m_L(2, p) * (m * val2 * (m_alpha[0][1] - m_alpha[1][1])),
						val * (m_L(1, p) * m_L(2, p) * m_alpha[0][2] + m_L(0, p) * m_L(2, p) * m_alpha[1][2] + m_L(0, p) * m_L(1, p) * m_alpha[2][2]) +
						m_L(0, p) * m_L(1, p) * m_L(2, p) * (m * val2 * (m_alpha[0][2] - m_alpha[1][2])));
	}
	val2 = std::pow(m_L(0, p) - m_L(1, p), m - 1);
	valp2 = std::pow(2 * m_L(2, p) - 1, n - 1);
	return Point(	val * valp * (m_L(1, p) * m_L(2, p) * m_alpha[0][1] + m_L(0, p) * m_L(2, p) * m_alpha[1][1] + m_L(0, p) * m_L(1, p) * m_alpha[2][1]) +
					m_L(0, p) * m_L(1, p) * m_L(2, p) * (m * val2 * valp * (m_alpha[0][1] - m_alpha[1][1]) + 2 * n * val * valp2 * m_alpha[2][1]), 
					val * valp * (m_L(1, p) * m_L(2, p) * m_alpha[0][2] + m_L(0, p) * m_L(2, p) * m_alpha[1][2] + m_L(0, p) * m_L(1, p) * m_alpha[2][2]) +
					m_L(0, p) * m_L(1, p) * m_L(2, p) * (m * val2 * valp * (m_alpha[0][2] - m_alpha[1][2]) + 2 * n * val * valp2 * m_alpha[2][2]));
}
//const function<const double(const Point&)> CTriangleLinearBasis::GetShapeFunction(const int k) const
//{
	//return [&](const Point& p){return m_alpha[k][0] + m_alpha[k][1] * p.x + m_alpha[k][2] * p.y;};
//}

CTriangleBasis::CTriangleBasis(const CTriangleBasis& t)
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
	m_number = t.m_number;
	m_order = t.m_order;
	m_1dorder = t.m_1dorder;
	//m_w = t.m_w;
	m_s = t.m_s;
	m_sp = t.m_sp;
	m_all = t.m_all;
}

void CTriangleBasis::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
{
	double pz;
	//m_det = (p2.y - p1.y)*(p3.z - p1.z) - (p3.y - p1.y)*(p2.z - p1.z);
	//m_det *= m_det;
	//py = (p3.x - p1.x)*(p2.z - p1.z) - (p2.x - p1.x)*(p3.z - p1.z);
	//py *= py;
	pz = (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
	//pz *= pz;
	//m_det = sqrt(m_det + py + pz);
	m_det = pz;
}

void CTriangleBasis::compAlpha(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
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

void CTriangleBasis::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
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

const double CTriangleBasis::GetValue(const Point& p) const
{
	return 0.;
}
const int CTriangleBasis::IncreaseOrder()
{
	++m_order;
	++m_1dorder;
	createS();
	m_all.push_back(m_sp);
	return 0;
}

const Point mid_point(const Point& p1, const Point& p2)
{
	return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
}

const Point s_point(const Point& p1, const Point& p2, const double s)
{
	return Point((p1.x + p2.x) / s, (p1.y + p2.y) / s);
}

const Point center_point(const Point& p1, const Point& p2, const Point& p3)
{
	return Point((p1.x + p2.x + p3.x) / 3, (p1.y + p2.y + p3.y) / 3);
}

const double CTriangleBasis::GetWeight(const int node, const vector<Point>& pts, const std::function<const double(const Point&)>& func) const
{
    //if (node < 3)
		//if(node < pts.size())
            //return func(pts[node]);

	vector<double> qq(15);
	qq[0] = func(pts[0]);
	qq[1] = func(pts[1]);
	qq[2] = func(pts[2]);
    qq[3] = 4 * func(mid_point(pts[0], pts[1])) - 2 * (func(pts[0]) + func(pts[1]));
    qq[4] = 4 * func(mid_point(pts[1], pts[2])) - 2 * (func(pts[1]) + func(pts[2]));
    qq[5] = 4 * func(mid_point(pts[0], pts[2])) - 2 * (func(pts[0]) + func(pts[2]));
	const Point center = Point{ (pts[0].x + pts[1].x + pts[2].x) / 3, (pts[0].y + pts[1].y + pts[2].y) / 3 };
	//const double q6 = 27. / 2 * func(Point{ pts[0].x + (pts[1].x - pts[0].x) / 3, pts[0].y + (pts[1].y - pts[0].y) / 3 }) - 9. * func(pts[0]) - 9. / 2 * func(pts[1]) - 3 * q3;
    qq[6] = 27. / 2 * (func(Point{pts[0].x  + (pts[1].x - pts[0].x) / 3, pts[0].y + (pts[1].y - pts[0].y) / 3}) - 2. / 3 * func(pts[0]) - 1. / 3 * func(pts[1]) - 2. / 9 * qq[3]);
    qq[7] = 27. / 2 * (func(Point{pts[1].x  + (pts[2].x - pts[1].x) / 3, pts[1].y + (pts[2].y - pts[1].y) / 3}) - 2. / 3 * func(pts[1]) - 1. / 3 * func(pts[2]) - 2. / 9 * qq[4]);
    qq[8] = 27. / 2 * (func(Point{pts[0].x  + (pts[2].x - pts[0].x) / 3, pts[0].y + (pts[2].y - pts[0].y) / 3}) - 2. / 3 * func(pts[0]) - 1. / 3 * func(pts[2]) - 2. / 9 * qq[5]);
    /*qq[6] = 32. / 3. * (func(Point{ pts[0].x + (pts[1].x - pts[0].x) / 4., pts[0].y + (pts[1].y - pts[0].y) / 4. }) - 0.75 * func(pts[0]) - 0.25 * func(pts[1]) - 3. / 16. * qq[3]);
	qq[7] = 32. / 3. * (func(Point{ pts[1].x + (pts[2].x - pts[1].x) / 4, pts[1].y + (pts[2].y - pts[1].y) / 4 }) - 0.75 * func(pts[1]) - 0.25 * func(pts[2]) - 3. / 16. * qq[4]);
    qq[8] = 32. / 3. * (func(Point{ pts[0].x + (pts[2].x - pts[0].x) / 4, pts[0].y + (pts[2].y - pts[0].y) / 4 }) - 0.75 * func(pts[0]) - 0.25 * func(pts[2]) - 3. / 16. * qq[5]);*/
	//const double q7 = 27. / 2 * func(Point{ pts[1].x + (pts[2].x - pts[1].x) / 3, pts[1].y + (pts[2].y - pts[1].y) / 3 }) - 9. * func(pts[1]) - 9. / 2 * func(pts[2]) - 3 * q4;
	//const double q8 = 27. / 2 * func(Point{ pts[0].x + (pts[2].x - pts[0].x) / 3, pts[0].y + (pts[2].y - pts[0].y) / 3 }) - 9. * func(pts[0]) - 9. / 2 * func(pts[2]) - 3 * q5;

    qq[10] = 64./3*(func(Point{ pts[0].x + (pts[1].x - pts[0].x) / 4, pts[0].y + (pts[1].y - pts[0].y) / 4 }) - 0.75 * func(pts[0]) - 0.25 * func(pts[1]) -
		3. / 16 * qq[3] - 3. / 32 * qq[6]);
    qq[11] = 64./3*(func(Point{ pts[1].x + (pts[2].x - pts[1].x) / 4, pts[1].y + (pts[2].y - pts[1].y) / 4 }) - 0.75 * func(pts[1]) - 0.25 * func(pts[2]) -
		3. / 16 * qq[4] - 3. / 32 * qq[7]);
    qq[12] = 64./3*(func(Point{ pts[0].x + (pts[2].x - pts[0].x) / 4, pts[0].y + (pts[2].y - pts[0].y) / 4 }) - 0.75 * func(pts[0]) - 0.25 * func(pts[2]) -
		3. / 16 * qq[5] - 3. / 32 * qq[8]);

    //qq[13] = 10;
    //qq[14] = 10;
    Point temp3{pts[2].x + (pts[1].x - pts[2].x) / 4, pts[2].y + (pts[1].y - pts[2].y) / 4};
    Point temp4{pts[0].x + (pts[1].x - pts[0].x) / 4, pts[0].y + (pts[1].y - pts[0].y) / 4};

    Point temp1{mid_point(pts[1], pts[2])};
    Point temp2{mid_point(pts[0], pts[2])};

    Point inter;
    inter.x = (temp1.x * temp2.y - temp1.y * temp2.x) * (temp3.x - temp4.x) - (temp1.x - temp2.x) * (temp3.x * temp4.y - temp3.y * temp4.x);
    inter.x /= (temp1.x - temp2.x) * (temp3.y - temp4.y) - (temp1.y - temp2.y) * (temp3.x - temp4.x);

    inter.y = (temp1.x * temp2.y - temp1.y * temp2.x) * (temp3.y - temp4.y) - (temp1.y - temp2.y) * (temp3.x * temp4.y - temp3.y * temp4.x);
    inter.y /= (temp1.x - temp2.x) * (temp3.y - temp4.y) - (temp1.y - temp2.y) * (temp3.x - temp4.x);
    double sum = 0.;
    for (auto i = 0; i < 13; ++i)
        sum += qq[i] * GetShapeFunction(i, inter);
    qq[9] = (func(inter) - sum) / GetShapeFunction(9, inter);
    //cout << GetShapeFunction(0, inter) << "\t" << GetShapeFunction(1, inter) << "\t" << GetShapeFunction(2, inter) << endl;
    temp3 = {pts[2].x + (pts[1].x - pts[2].x) / 8, pts[2].y + (pts[1].y - pts[2].y) / 8};
    temp4 = {pts[0].x + (pts[1].x - pts[0].x) / 8, pts[0].y + (pts[1].y - pts[0].y) / 8};
    inter.x = (temp1.x * temp2.y - temp1.y * temp2.x) * (temp3.x - temp4.x) - (temp1.x - temp2.x) * (temp3.x * temp4.y - temp3.y * temp4.x);
    inter.x /= (temp1.x - temp2.x) * (temp3.y - temp4.y) - (temp1.y - temp2.y) * (temp3.x - temp4.x);

    inter.y = (temp1.x * temp2.y - temp1.y * temp2.x) * (temp3.y - temp4.y) - (temp1.y - temp2.y) * (temp3.x * temp4.y - temp3.y * temp4.x);
    inter.y /= (temp1.x - temp2.x) * (temp3.y - temp4.y) - (temp1.y - temp2.y) * (temp3.x - temp4.x);


    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0., 1.);
    double r1, r2;
    r1 = dist(mt);
    r2 = dist(mt);

    //cout << GetShapeFunction(0, inter) << "\t" << GetShapeFunction(1, inter) << "\t" << GetShapeFunction(13, inter) << endl;
    sum = 0.;
    for (auto i = 0; i < 13; ++i)
        sum += qq[i] * GetShapeFunction(i, inter);
    qq[14] = (func(inter) - sum) / GetShapeFunction(14, inter);
    sum = 0;
    for (auto i = 0; i < 13; ++i)
        sum += qq[i] * GetShapeFunction(i, center);
    qq[13] = (func(center) - sum) / GetShapeFunction(13, center);
    /*Point tp{inter};
    tp = Point((1 - sqrt(r1))*pts[0].x + sqrt(r1)*(1 - r2)*pts[1].x + sqrt(r1)*r2*pts[2].x, (1 - sqrt(r1))*pts[0].y + sqrt(r1)*(1 - r2)*pts[1].y + sqrt(r1)*r2*pts[2].y);
    cout << GetShapeFunction(0, tp) << "\t" << GetShapeFunction(1, tp) << "\t" << GetShapeFunction(2, tp) << endl;
    //const Point tp(Point{ pts[0].x + (pts[1].x - pts[0].x) / 1, pts[0].y + (pts[1].y - pts[0].y) / 1 });
    double expected = func(tp);
    Point dexpected(4. * tp.x * tp.x * tp.x, 4. * tp.y * tp.y * tp.y);
    double actual = 0;
    Point dactual;
    for (auto i = 0; i < m_s; ++i)
        actual += qq[i] * GetShapeFunction(i, tp);

    for (auto i = 0; i < m_s; ++i)
    {
        dactual.x += qq[i] * GetGradShapeFunction(i, tp).x;
        dactual.y += qq[i] * GetGradShapeFunction(i, tp).y;
    }
    cout << node <<  "\t" << actual << "\t" << expected << endl;
    cout << dactual.x << "\t" << dexpected.x << endl;
    cout << dactual.y << "\t" << dexpected.y << endl;*/
    //if(m_s < 10)
        return qq[node];
    //return 0.;
	const Point c01 = Point{ (pts[0].x + pts[1].x) / 2, (pts[0].y + pts[1].y) / 2 };
	const Point c02 = Point{ (pts[0].x + pts[2].x) / 2, (pts[0].y + pts[2].y) / 2 };
	const Point c12 = Point{ (pts[1].x + pts[2].x) / 2, (pts[1].y + pts[2].y) / 2 };


	Matrix matrix;
	const double m_ss = m_s - 6;
	matrix.Create(m_s);
	vector<double> qq2(m_s);
	Point temp;


	matrix.AddElement(0, 0, 1);
	matrix.AddElement(1, 1, 1);
	matrix.AddElement(2, 2, 1);
	matrix.AddElement(3, 3, 1);
	matrix.AddElement(4, 4, 1);
	matrix.AddElement(5, 5, 1);
	qq2[0] = qq[0];
	qq2[1] = qq[1];
	qq2[2] = qq[2];
	qq2[3] = qq[3];
	qq2[4] = qq[4];
	qq2[5] = qq[5];
	for (auto i = 6; i < m_s; ++i)
	{
		r1 = dist(mt);
		r2 = dist(mt);
		temp = Point((1 - sqrt(r1))*pts[0].x + sqrt(r1)*(1 - r2)*pts[1].x + sqrt(r1)*r2*pts[2].x, (1 - sqrt(r1))*pts[0].y + sqrt(r1)*(1 - r2)*pts[1].y + sqrt(r1)*r2*pts[2].y);
		for (auto j = 0; j < m_s; ++j)
			matrix.AddElement(i, j, GetShapeFunction(j, temp));
		qq2[i] = func(temp);
	}
    //return 0.;
    ESolver solver;
	vector<double> sol(m_s);
	solver.Gauss(matrix, qq2, sol);
    /*expected = func(tp);
    actual = 0;
    Point pact;
	for (auto i = 0; i < m_s; ++i)
        actual += sol[i] * GetShapeFunction(i, tp);
    cout << node <<  "\t" << actual << "\t" << expected << endl;

    for (auto i = 0; i < m_s; ++i)
    {
        pact.x += sol[i] * GetGradShapeFunction(i, tp).x;
        pact.y += sol[i] * GetGradShapeFunction(i, tp).y;
    }
    cout << pact.x << "\t" << dexpected.x << endl;
    cout << pact.y << "\t" << dexpected.y << endl;
    cout << "----------------" << endl;
    for (auto i = 0; i < m_s; ++i)
    {
        cout << i << "\t" << qq[i] << "\t" << sol[i] << endl;
    }
    cout << "----------------" << endl;*/
    return sol[node];
}
//const int CTriangleBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const CFESolution CTriangleBasis::GetValue(const int n) const
//{
//	return m_w[n];
//}








































