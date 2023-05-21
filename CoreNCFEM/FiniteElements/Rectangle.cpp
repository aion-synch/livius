#include "Rectangle.h"
using namespace std;
using namespace corenc;
using namespace Mesh;

CRectangle::CRectangle()
{
	m_edges[0] = m_edges[1] = m_edges[2] = -1;
	m_nodes.resize(4);
	m_order = 0;
	m_number = 0;
    m_px = m_py = 0;
}

CRectangle::CRectangle(const int n1, const int n2, const int n3, const int n4, const int order)
{
	m_nodes.resize(4);
	m_nodes[0] = n1;
	m_nodes[1] = n2;
	m_nodes[2] = n3;
	m_nodes[3] = n4;
	m_order = order;
	m_number = 1;
	m_edges[0] = m_edges[1] = m_edges[2] = m_edges[3] = -1;
    m_px = m_py = 1;
}

CRectangle::CRectangle(const int n1, const int n2, const int n3, const int n4, const int e1, const int e2, const int e3, const int e4, const int order)
{
	m_nodes.resize(4);
	m_nodes[0] = n1;
	m_nodes[1] = n2;
	m_nodes[2] = n3;
	m_nodes[3] = n4;
	m_edges[0] = e1;
	m_edges[1] = e2;
	m_edges[2] = e3;
	m_edges[3] = e4;
	m_order = order;
    m_px = m_py = 1;
	m_number = 1;
}

CRectangle::CRectangle(const int* nodes, const int order)
{
	m_order = order;
	m_number = 1;
	m_nodes.resize(order);
	for (int i = 0; i < m_order; ++i)
		m_nodes[i] = nodes[i];
    m_px = m_py = 1;
}

CRectangle::CRectangle(const int* nodes, const int* edges, const int order)
{
	m_order = order;
	m_nodes.resize(order);
	for (int i = 0; i < m_order; ++i)
		m_nodes[i] = nodes[i];
	m_edges[0] = edges[0];
	m_edges[1] = edges[1];
	m_edges[2] = edges[2];
	m_edges[3] = edges[3];
	m_number = 1;
    m_px = m_py = 1;
}

CRectangle::CRectangle(const CRectangle& t)
{
	vector<int>().swap(m_nodes);
	m_nodes = t.m_nodes;
	m_order = t.m_order;
	m_edges[0] = t.m_edges[0];
	m_edges[1] = t.m_edges[1];
	m_edges[2] = t.m_edges[2];
	m_edges[3] = t.m_edges[3];
	m_number = t.m_number;
    m_px = t.m_px;
    m_py = t.m_py;
}

const int CRectangle::SetOrder(const int px, const int py)
{
    m_px = px;
    m_py = py;
    m_order = (px + 1) * (py + 1);
    m_nodes.resize(m_order);
    return 0;
}

const int CRectangle::GetNode(const int n) const
{
	if (n < m_nodes.size())
		return m_nodes[n];
	cout << "GetNode(" << n << "): Wrong number of the node." << endl;
	return -1;
}

const int CRectangle::GetNode(const NODES& node) const
{
	if (node == NODES::FIRST)
		return m_nodes[0];
	return m_nodes[1];
}

const int CRectangle::GetEdge(const int n) const
{
	if (n < 4)
		return m_edges[n];
	cout << "GetEdge(" << n << "): Wrong number of the edge." << endl;
	return -1;
}

const int CRectangle::GetFacet(const int) const
{
	return -1;
}

const int CRectangle::GetNumberOfNodes() const
{
	return (int)m_nodes.size();
}

const int CRectangle::GetNumberOfEdges() const
{
	return 4;
}

const int CRectangle::GetNumberOfFacets() const
{
	return 1;
}

void CRectangle::SetNode(const int n, const int node)
{
	if (n < m_order)
		m_nodes[n] = node;
	else
		cout << "SetNode(" << n << ", " << node << "): Wrong number of the node." << endl;
}

const int CRectangle::IncreaseOrder()
{
	++m_number;
	//const auto pr = m_order;
	m_order = 3 * (m_number)+(m_number - 1) * (m_number - 2) / 2;
	m_nodes.resize(m_order);
	//for (int i = pr; i < m_order; ++i)
	//m_nodes[i] = -1;
	return 0;
}

void CRectangle::SetEdge(const int k, const int edge)
{
	if (k < 4)
		m_edges[k] = edge;
	else
		cout << "SetEdge(" << k << ", " << edge << "): Wrong number of the edge." << endl;
}

void CRectangle::SetFacet(const int, const int)
{
	cout << "You tried." << endl;
}

const double CRectangle::Integrate(const std::function<const double(const Point &)> &f, const std::vector<Point>& v) const
{
	Point temp;
	double integral = 0;
	double ax{ (v[3].x - v[0].x) / 2 };
	double ay{ (v[3].y - v[0].y) / 2 };
	if (fabs(ax) < 1e-11)
	{
		ax = (v[3].z - v[0].z) / 2;
		for (int i = 0; i < 12; ++i)
		{
			temp.x = v[0].x;
			temp.z = v[0].z + ax * (GaussRectangular::m_ra[i] + 1);
			temp.y = v[0].y + ay * (GaussRectangular::m_rb[i] + 1);
			integral += GaussRectangular::m_rw[i] * f(temp)*temp.Jacobian();
		}
		return integral;
	}
	if (fabs(ay) < 1e-11)
	{
		ay = (v[3].z - v[0].z) / 2;
		for (int i = 0; i < 12; ++i)
		{
			temp.x = v[0].x + ax * (GaussRectangular::m_ra[i] + 1);
			temp.z = v[0].z + ay * (GaussRectangular::m_rb[i] + 1);
			temp.y = v[0].y;
			integral += GaussRectangular::m_rw[i] * f(temp)*temp.Jacobian();
		}
		return integral;
	}
	for (int i = 0; i < 12; ++i)
	{
		temp.x = v[0].x + ax * (GaussRectangular::m_ra[i] + 1);
		temp.y = v[0].y + ay * (GaussRectangular::m_rb[i] + 1);
		temp.z = v[0].z;
		integral += GaussRectangular::m_rw[i] * f(temp)*temp.Jacobian();
	}
	return ax * ay * integral;
}

const Point CRectangle::Integrate(const std::function<const Point(const Point &)> &f, const std::vector<Point>& v) const
{
	// triangle!!
	Point temp;
	Point integral{ 0,0,0 };
	for (int i = 0; i < 16; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		integral += GaussTriangle::m_trw[i] * f(temp)*temp.Jacobian();
	}
	return integral;
}

const vector<double> CRectangle::Integrate(const std::function<const std::vector<double>(const Point &)> &f, const std::vector<Point> &v) const
{
	// triangle !!
	Point temp;
	auto size{ f({ 0,0,0 }).size() };
	vector<double> integral;
	integral.resize(size);
	for (int i = 0; i < 16; ++i)
	{
		temp.x = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].x + GaussTriangle::m_tra[i] * v[1].x + GaussTriangle::m_trb[i] * v[2].x;
		temp.y = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].y + GaussTriangle::m_tra[i] * v[1].y + GaussTriangle::m_trb[i] * v[2].y;
		temp.z = (1 - GaussTriangle::m_tra[i] - GaussTriangle::m_trb[i]) * v[0].z + GaussTriangle::m_tra[i] * v[1].z + GaussTriangle::m_trb[i] * v[2].z;
		auto fx{ f(temp) };
		for (int j = 0; j < size; ++j)
			integral[j] = GaussTriangle::m_trw[i] * fx[j] * temp.Jacobian();
	}
	return integral;
}

CRectangleBasis::CRectangleBasis()
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
CRectangleBasis::CRectangleBasis(const Point& p1, const Point& p2, const Point& p3, const Point& p4, const int order)
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
	m_number = 4;
}

CRectangleBasis::CRectangleBasis(const Point* p, const int order)
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
}

const int CRectangleBasis::createS()
{
	m_sp = 4;
	m_s = 4;
	return 0;
}

const int CRectangleBasis::GetNumberOfShapeFunctions() const
{
    return m_sp;
}

const Point CRectangleBasis::GetNormal() const
{
	return m_normal;
}

void CRectangleBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}

const double CRectangleBasis::m_x1(const double x) const
{
	return (m_points[3].x - x) / m_hx;
}

const double CRectangleBasis::m_x2(const double x) const
{
	return (x - m_points[0].x) / m_hx;
}

const double CRectangleBasis::m_y1(const double y) const
{
	return (m_points[3].y - y) / m_hy;
}

const double CRectangleBasis::m_y2(const double y) const
{
	return (y - m_points[0].y) / m_hy;
}

const double CRectangleBasis::GetShapeFunction(const int k, const Point &p) const
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
	return 0;
}

const Point CRectangleBasis::GetGradShapeFunction(const int k, const Point & p) const
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
	default:
		break;
	}
	return m_dpsi;
}

CRectangleBasis::CRectangleBasis(const CRectangleBasis& t)
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

void CRectangleBasis::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
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

void CRectangleBasis::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
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

const double CRectangleBasis::GetValue(const Point& p) const
{
	return 0.;
}
const int CRectangleBasis::IncreaseOrder()
{
	createS();
	++m_order;
	++m_1dorder;
	//m_all.push_back(m_sp);
	return 0;
}
const double CRectangleBasis::GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
{
    //return 0.0;
    return f(verts[node]);
}
//const int CRectangleBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const CFESolution CRectangleBasis::GetValue(const int n) const
//{
//	return m_w[n];
//}


CRectangleConstantBasis::CRectangleConstantBasis()
{
	m_det = 0;
	createS();
	m_1dorder = m_order = 0;
	m_hx = 1;
	m_hy = 1;
	m_points.resize(4);
	m_normal = Point(0, 0, 0);
}
CRectangleConstantBasis::CRectangleConstantBasis(const Point& p1, const Point& p2, const Point& p3, const Point& p4, const int order)
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
	m_number = 4;
}

CRectangleConstantBasis::CRectangleConstantBasis(const Point* p, const int order)
{
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
}

const int CRectangleConstantBasis::createS()
{
	m_sp = 1;
	m_s = 1;
	return 0;
}

const int CRectangleConstantBasis::GetNumberOfShapeFunctions() const
{
	return m_sp;
}

const Point CRectangleConstantBasis::GetNormal() const
{
	return m_normal;
}

void CRectangleConstantBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}

const double CRectangleConstantBasis::m_x1(const double x) const
{
	return (m_points[3].x - x) / m_hx;
}

const double CRectangleConstantBasis::m_x2(const double x) const
{
	return (x - m_points[0].x) / m_hx;
}

const double CRectangleConstantBasis::m_y1(const double y) const
{
	return (m_points[3].y - y) / m_hy;
}

const double CRectangleConstantBasis::m_y2(const double y) const
{
	return (y - m_points[0].y) / m_hy;
}

const double CRectangleConstantBasis::GetShapeFunction(const int k, const Point &p) const
{
	return 1.;
}

const Point CRectangleConstantBasis::GetGradShapeFunction(const int k, const Point & p) const
{
	return Point(0, 0, 0);
}

CRectangleConstantBasis::CRectangleConstantBasis(const CRectangleConstantBasis& t)
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

void CRectangleConstantBasis::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
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

void CRectangleConstantBasis::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
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

const double CRectangleConstantBasis::GetValue(const Point& p) const
{
	return 0.;
}
const int CRectangleConstantBasis::IncreaseOrder()
{
	createS();
	++m_order;
	++m_1dorder;
	//m_all.push_back(m_sp);
	return 0;
}
//const int CRectangleConstantBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const CFESolution CRectangleConstantBasis::GetValue(const int n) const
//{
//	return m_w[n];
//}
