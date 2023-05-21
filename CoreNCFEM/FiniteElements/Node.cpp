#include "Node.h"
using namespace corenc;
using namespace Mesh;
using namespace std;

CNode::CNode()
{
	m_node = -1;
}

CNode::CNode(const CNode& n)
{
	m_node = n.m_node;
}

CNode::CNode(const int* n)
{
	m_node = n[0];
}

const int CNode::GetNode(const int n) const
{
	if(n < 1)
		return m_node;
	return -1;
}

const int CNode::GetNode(const NODES& node) const
{
	return m_node;
}
const int CNode::GetNumberOfNodes() const
{
	return 1;
}

void CNode::SetNode(const int, const int node)
{
	m_node = node;
}

const double CNode::Integrate(const std::function<const double (const Point &)> &f, const std::vector<Point> &v) const
{
	return f(v[0]);
}

const vector<double> CNode::Integrate(const std::function<const vector<double> (const Point &)> &f, const std::vector<Point> &v) const
{
	return f(v[0]);
}

const Point CNode::Integrate(const std::function<const Point (const Point &)> &f, const std::vector<Point> &v) const
{
	return f(v[0]);
}

CNodeBasis::CNodeBasis()
{
	m_p0 = Point{0,0,0};
	m_normal = Point{1,0,0};
}


CNodeBasis::CNodeBasis(const Point* p)
{
	m_p0 = p[0];
	m_normal = Point{1,0,0};
}

const int CNodeBasis::GetNumberOfShapeFunctions() const
{
	return 1;
}

const double CNodeBasis::GetShapeFunction(const int, const Mesh::Point &) const
{
	return m_p0.x;
}

const Point CNodeBasis::GetGradShapeFunction(const int, const Mesh::Point &) const
{
	return Point{0,0,0};
}

const Point CNodeBasis::GetNormal() const
{
	return m_normal;
}

void CNodeBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}
//const CFESolution CNodeBasis::GetValue(const Point& p) const
//{
//	return m_w;
//}
//const int CNodeBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const int CNodeBasis::SetValue(const int, const CFESolution& value)
//{
//	return 0;
//}
//const CFESolution CNodeBasis::GetValue(const int) const
//{
//	return m_w;
//}


