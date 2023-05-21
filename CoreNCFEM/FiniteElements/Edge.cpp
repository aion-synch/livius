#include "Edge.h"
#include <iostream>
using namespace corenc;
using namespace Mesh;
using namespace std;

CEdge::CEdge()
{
	m_number = 2;
	m_nodes.resize(2);
	m_nodes[0] = m_nodes[1] = -1;
}

CEdge::CEdge(const CEdge& e)
{
	vector<int>().swap(m_nodes);
	m_nodes = e.m_nodes;
	m_number = e.m_number;
}

CEdge::CEdge(const int n1, const int n2)
{
	m_number = 2;
	m_nodes.resize(2);
	m_nodes[0] = n1;
	m_nodes[1] = n2;
}

CEdge::CEdge(const int* n)
{
	m_number = 2;
	m_nodes.resize(2);
	m_nodes[0] = n[0];
	m_nodes[1] = n[1];
}

const double CEdge::Integrate(const std::function<const double (const Point &)> &f, const std::vector<Point>&p) const
{
	double m_mes = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
	double integral = 0;
	Point temp;
	double param{ sqrt(((p[1].x - p[0].x)*(p[1].x - p[0].x) / (m_mes*m_mes) + (p[1].y - p[0].y)*(p[1].y - p[0].y) / (m_mes*m_mes))) };
	for (int k = 0; k < Gauss1dim::m_order; ++k)
	{
		//temp.x = p[0].x + 0.5*(m_mes + m_mes * Gauss1dim::m_a[k]) * (p[1].x - p[0].x) / m_mes;
		temp.x = m_mes / 2 * Gauss1dim::m_a[k] + (p[0].x + p[1].x) / 2;
		temp.y = p[0].y + 0.5*(m_mes + m_mes * Gauss1dim::m_a[k]) * (p[1].y - p[0].y) / m_mes;
		integral += Gauss1dim::m_w[k] * f(temp)*temp.Jacobian();
	}
	//integral *= m_mes*param / 2;// sqrt((p[1].x - p[0].x) / h + (p[1].y - p[0].y) / h);
	integral *= m_mes / 2;
	return integral;
}

const vector<double> CEdge::Integrate(const std::function<const vector<double> (const Point &)> &f, const std::vector<Point>&p) const
{
	double m_mes = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
	auto size {f({0,0,0}).size()};
	vector<double> integral;
	integral.resize(size);
	Point temp;
	double param{ sqrt(((p[1].x - p[0].x)*(p[1].x - p[0].x) / (m_mes*m_mes) + (p[1].y - p[0].y)*(p[1].y - p[0].y) / (m_mes*m_mes))) };
	for (int k = 0; k < Gauss1dim::m_order; ++k)
	{
		temp.x = p[0].x + 0.5*(m_mes + m_mes * Gauss1dim::m_a[k]) * (p[1].x - p[0].x) / m_mes;
		temp.y = p[0].y + 0.5*(m_mes + m_mes * Gauss1dim::m_a[k]) * (p[1].y - p[0].y) / m_mes;
		auto fx{f(temp)};
		for(int j = 0; j < size; ++j)
			integral[j] += Gauss1dim::m_w[k] * fx[j]*temp.Jacobian();
	}
	for(int j = 0; j < size; ++j)
		integral[j] *= m_mes*param / 2;// sqrt((p[1].x - p[0].x) / h + (p[1].y - p[0].y) / h);
	return integral;
}

const Point CEdge::Integrate(const std::function<const Point (const Point &)> &f, const std::vector<Point>&p) const
{
	double m_mes = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
	Point integral{0,0,0};
	Point temp;
	double param{ sqrt(((p[1].x - p[0].x)*(p[1].x - p[0].x) / (m_mes*m_mes) + (p[1].y - p[0].y)*(p[1].y - p[0].y) / (m_mes*m_mes))) };
	for (int k = 0; k < Gauss1dim::m_order; ++k)
	{
		temp.x = p[0].x + 0.5*(m_mes + m_mes * Gauss1dim::m_a[k]) * (p[1].x - p[0].x) / m_mes;
		temp.y = p[0].y + 0.5*(m_mes + m_mes * Gauss1dim::m_a[k]) * (p[1].y - p[0].y) / m_mes;
		integral += Gauss1dim::m_w[k] * f(temp)*temp.Jacobian();
	}
	integral *= m_mes*param / 2;// sqrt((p[1].x - p[0].x) / h + (p[1].y - p[0].y) / h);
	return integral;
}

const int CEdge::GetNode(const int k) const
{
	if(k < m_nodes.size())
		return m_nodes[k];
	cout << "CEdge::GetNode(" << k << "): k is out of range." << endl;
	return -1;
}

const int CEdge::GetNode(const NODES& node) const
{
	if (node == NODES::FIRST)
		return m_nodes[0];
	return m_nodes[1];
}
const int CEdge::IncreaseOrder()
{
	++m_number;
	m_nodes.resize(m_number);
	return 0;
}
const int CEdge::GetNumberOfNodes() const
{
	return m_number;
}

void CEdge::SetNode(const int k, const int node)
{
	if (k < m_nodes.size())
	{
		m_nodes[k] = node;
		return;
	}
	cout << "CEdge::SetNode(" << k << ", " << node << "): k is out of range." << endl;
}

CEdgeLinearBasis::CEdgeLinearBasis()
{
	m_number = 2;
	m_p0 = Point(0, 0, 0);
	m_p1 = Point(0, 0, 0);
	m_normal = Point(0, 0, 0);
	m_mes = 0;
}

CEdgeLinearBasis::CEdgeLinearBasis(const Point& p0, const Point& p1)
{
	m_number = 2;
	m_p0 = p0;
	m_p1 = p1;
	CompLenght();
	CompNormal();
}

CEdgeLinearBasis::CEdgeLinearBasis(const Point* p)
{
	m_number = 2;
	m_p0 = p[0];
	m_p1 = p[1];
	CompLenght();
	CompNormal();
}
CEdgeLinearBasis::CEdgeLinearBasis(const CEdgeLinearBasis& e)
{
	m_number = e.m_number;
	m_p0 = e.m_p0;
	m_p1 = e.m_p1;
	m_normal = e.m_normal;
	m_mes = e.m_mes;
}

void CEdgeLinearBasis::CompNormal()
{
	m_normal.x = m_p1.y - m_p0.y;
	m_normal.y = m_p1.x - m_p0.x;
	const double m = sqrt(m_normal.x * m_normal.x + m_normal.y * m_normal.y);
	m_normal.x /= m;
	m_normal.y /= m;
}
const int CEdgeLinearBasis::IncreaseOrder() 
{ 
	++m_number; 
	return 0;
};
void CEdgeLinearBasis::CompLenght()
{
	m_mes = sqrt((m_p1.x - m_p0.x)*(m_p1.x - m_p0.x) + (m_p1.y - m_p0.y)*(m_p1.y - m_p0.y) + (m_p1.z - m_p0.z)*(m_p1.z - m_p0.z));
}

const double CEdgeLinearBasis::GetShapeFunction(const int k, const Mesh::Point &p) const
{
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes};
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	if(p.x < m_p0.x || p.x > m_p1.x)
		return 0;
	switch (k)
	{
		case 0:
			return 1 - temp.x;
			//return 1;
		case 1:
			//return 2*temp.x-1;
			return temp.x;
		default:
			return -1;
	}
}

//const CFESolution CEdgeLinearBasis::GetValue(const Point& p) const
//{
//	return m_w[0];
//}
//const int CEdgeLinearBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const int CEdgeLinearBasis::SetValue(const int, const CFESolution& value)
//{
//	return 0;
//}
//const CFESolution CEdgeLinearBasis::GetValue(const int n) const
//{
	//return dynamic_cast<CSolution*>(&m_w[0]);
//	return m_w[n];
//}

const Point CEdgeLinearBasis::GetGradShapeFunction(const int k, const Mesh::Point&) const
{
	if(k)
		return Point(1/m_mes, 0, 0);
	return Point(-1/m_mes, 0, 0);
	return Point(0,0,0);
}

const int CEdgeLinearBasis::GetNumberOfShapeFunctions() const
{
	return m_number;
}

const Point CEdgeLinearBasis::GetNormal() const
{
	return m_normal;
}

void CEdgeLinearBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}



CEdgeConstantBasis::CEdgeConstantBasis()
{
	m_p0 = Point(0, 0, 0);
	m_p1 = Point(0, 0, 0);
	m_normal = Point(0, 0, 0);
	m_mes = 0;
}

CEdgeConstantBasis::CEdgeConstantBasis(const Point& p0, const Point& p1)
{
	m_p0 = p0;
	m_p1 = p1;
	CompLenght();
	CompNormal();
}

CEdgeConstantBasis::CEdgeConstantBasis(const Point* p)
{
	m_p0 = p[0];
	m_p1 = p[1];
	CompLenght();
	CompNormal();
}
CEdgeConstantBasis::CEdgeConstantBasis(const CEdgeConstantBasis& e)
{
	m_p0 = e.m_p0;
	m_p1 = e.m_p1;
	m_normal = e.m_normal;
	m_mes = e.m_mes;
}

void CEdgeConstantBasis::CompNormal()
{
	m_normal.x = m_p0.y - m_p1.y;
	m_normal.y = m_p1.x - m_p0.x;
}

void CEdgeConstantBasis::CompLenght()
{
	m_mes = sqrt((m_p1.x - m_p0.x)*(m_p1.x - m_p0.x) + (m_p1.y - m_p0.y)*(m_p1.y - m_p0.y) + (m_p1.z - m_p0.z)*(m_p1.z - m_p0.z));
}

const double CEdgeConstantBasis::GetShapeFunction(const int, const Mesh::Point &p) const
{
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes};
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	return 1;
	return m_mes;
}

const Point CEdgeConstantBasis::GetGradShapeFunction(const int, const Mesh::Point&) const
{
	return Point(0, 0, 0);
}

const int CEdgeConstantBasis::GetNumberOfShapeFunctions() const
{
	return m_number;
}

const Point CEdgeConstantBasis::GetNormal() const
{
	return m_normal;
}

void CEdgeConstantBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}

//const CFESolution CEdgeConstantBasis::GetValue(const Point& p) const
//{
//	double sum{ 0. };
//	return m_w;
//}
//const int CEdgeConstantBasis::SetValue(const int, CSolution* value)
//{
//	const CFESolution* a = static_cast<const CFESolution*>(value);
//	m_w = *a;
//	return 0;
//}
//const int CEdgeConstantBasis::SetValue(const int, const CFESolution& value)
//{
//	m_w = value;
//	return 0;
//}
//const CFESolution CEdgeConstantBasis::GetValue(const int n) const
//{
//	return m_w;
//}

CEdge2ndBasis::CEdge2ndBasis()
{
	m_p0 = Point(0, 0, 0);
	m_p1 = Point(0, 0, 0);
	m_normal = Point(0, 0, 0);
	m_mes = 0;
}

CEdge2ndBasis::CEdge2ndBasis(const Point& p0, const Point& p1)
{
	m_p0 = p0;
	m_p1 = p1;
	CompLenght();
	CompNormal();
}

CEdge2ndBasis::CEdge2ndBasis(const Point* p)
{
	m_p0 = p[0];
	m_p1 = p[1];
	CompLenght();
	CompNormal();
}
CEdge2ndBasis::CEdge2ndBasis(const CEdge2ndBasis& e)
{
	m_p0 = e.m_p0;
	m_p1 = e.m_p1;
	m_normal = e.m_normal;
	m_mes = e.m_mes;
}

void CEdge2ndBasis::CompNormal()
{
	m_normal.x = m_p0.y - m_p1.y;
	m_normal.y = m_p1.x - m_p0.x;
}

void CEdge2ndBasis::CompLenght()
{
	m_mes = sqrt((m_p1.x - m_p0.x)*(m_p1.x - m_p0.x) + (m_p1.y - m_p0.y)*(m_p1.y - m_p0.y) + (m_p1.z - m_p0.z)*(m_p1.z - m_p0.z));
}

const double CEdge2ndBasis::GetShapeFunction(const int k, const Mesh::Point &p) const
{
	if(p.x < m_p0.x || p.x > m_p1.x)
		return 0;
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes};
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	
	switch (k)
	{
		case 0:
			return 2 * (temp.x - 1) * (temp.x - 0.5);
		case 1:
			return -4 * temp.x * (temp.x - 1);
		case 2:
			return 2 * temp.x * (temp.x - 0.5);
		default:
			return -1;
	}
}

const Point CEdge2ndBasis::GetGradShapeFunction(const int k, const Mesh::Point &p) const
{
	if(p.x < m_p0.x || p.x > m_p1.x)
		return Point(0,0,0);
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes};
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	auto jac{1/m_mes};
	switch(k)
	{
		case 0:
			return Point(jac*(4 * temp.x - 3), 0, 0);
		case 1:
			return Point(jac*(-8 * temp.x + 4), 0, 0);
		case 2:
			return Point(jac*(4 * temp.x - 1), 0, 0);
		default:
			return Point(0,0,0);
	};
}

const int CEdge2ndBasis::GetNumberOfShapeFunctions() const
{
	return m_number;
}

const Point CEdge2ndBasis::GetNormal() const
{
	return m_normal;
}

void CEdge2ndBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}
//const CFESolution CEdge2ndBasis::GetValue(const Point& p) const
//{
//	return m_w[0];
//}
//const int CEdge2ndBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const int CEdge2ndBasis::SetValue(const int, const CFESolution& value)
//{
//	return 0;
//}
//const CFESolution CEdge2ndBasis::GetValue(const int n) const
//{
//	return m_w[n];
//}

CEdgeHermiteBasis::CEdgeHermiteBasis()
{
	//m_w.resize(4);
	m_p0 = Point(0, 0, 0);
	m_p1 = Point(0, 0, 0);
	m_normal = Point(0, 0, 0);
	m_mes = 0;
	
}

CEdgeHermiteBasis::CEdgeHermiteBasis(const Point& p0, const Point& p1)
{
	m_p0 = p0;
	m_p1 = p1;
	CompLenght();
	CompNormal();
	//m_w.resize(4);
}

CEdgeHermiteBasis::CEdgeHermiteBasis(const Point* p)
{
	m_p0 = p[0];
	m_p1 = p[1];
	CompLenght();
	CompNormal();
	//m_w.resize(4);
}
CEdgeHermiteBasis::CEdgeHermiteBasis(const CEdgeHermiteBasis& e)
{
	m_p0 = e.m_p0;
	m_p1 = e.m_p1;
	m_normal = e.m_normal;
	m_mes = e.m_mes;
	//m_w = e.m_w;
}

void CEdgeHermiteBasis::CompNormal()
{
	m_normal.x = m_p0.y - m_p1.y;
	m_normal.y = m_p1.x - m_p0.x;
}

void CEdgeHermiteBasis::CompLenght()
{
	m_mes = sqrt((m_p1.x - m_p0.x)*(m_p1.x - m_p0.x) + (m_p1.y - m_p0.y)*(m_p1.y - m_p0.y) + (m_p1.z - m_p0.z)*(m_p1.z - m_p0.z));
}

const double CEdgeHermiteBasis::GetShapeFunction(const int k, const Mesh::Point &p) const
{
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes };
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	if (p.x < m_p0.x || p.x > m_p1.x)
		return 0;
	switch (k)
	{
	case 0:
		return (1 - 3 * temp.x * temp.x + 2 * temp.x * temp.x * temp.x);
	case 1:
		return m_mes * (temp.x - 2 * temp.x * temp.x + temp.x * temp.x * temp.x);
	case 2:
		return 3 * temp.x * temp.x - 2 * temp.x * temp.x * temp.x;
	case 3:
		return m_mes * (temp.x * temp.x * temp.x - temp.x * temp.x);
	default:
		return -1;
	}
}

const Point CEdgeHermiteBasis::GetGradShapeFunction(const int k, const Mesh::Point& p) const
{
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes };
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	switch (k)
	{
	case 0:
		return Point((-6 * temp.x + 6 * temp.x * temp.x) / m_mes, 0, 0);
	case 1:
		return Point(1 - 4 * temp.x + 3 * temp.x * temp.x, 0, 0);
	case 2:
		return Point((6 * temp.x - 6 * temp.x * temp.x) / m_mes, 0, 0);
	case 3:
		return Point(-2 * temp.x + 3 * temp.x * temp.x, 0, 0);
	}
	return Point(0, 0, 0);
}

const int CEdgeHermiteBasis::GetNumberOfShapeFunctions() const
{
	return m_number;
}

const Point CEdgeHermiteBasis::GetNormal() const
{
	return m_normal;
}

void CEdgeHermiteBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}


//const CFESolution CEdgeHermiteBasis::GetValue(const Point& p) const
//{
//	return m_w[0] * GetShapeFunction(0, p) + m_w[1] * GetShapeFunction(1, p) + m_w[2] * GetShapeFunction(2, p) + m_w[3] * GetShapeFunction(3, p);
//}
//const int CEdgeHermiteBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const int CEdgeHermiteBasis::SetValue(const int n, const CFESolution& value)
//{
//	if (n < m_w.size())
//	{
//		m_w[n] = value;
//		return 0;
//	}
//	return 1;
//}
//const CFESolution CEdgeHermiteBasis::GetValue(const int n) const 
//{
//	if(n < m_w.size())
//		return m_w[n];
//	return CFESolution(0);
//}


CEdgeMultiBasis::CEdgeMultiBasis()
{
	//m_w.resize(2);
	m_p0 = Point(0, 0, 0);
	m_p1 = Point(0, 0, 0);
	m_normal = Point(0, 0, 0);
	m_mes = 0;
}

CEdgeMultiBasis::CEdgeMultiBasis(const Point& p0, const Point& p1)
{
	m_p0 = p0;
	m_p1 = p1;
	CompLenght();
	CompNormal();
	//m_w.resize(2);
}

CEdgeMultiBasis::CEdgeMultiBasis(const Point* p)
{
	m_p0 = p[0];
	m_p1 = p[1];
	CompLenght();
	CompNormal();
	//m_w.resize(2);
}
CEdgeMultiBasis::CEdgeMultiBasis(const CEdgeMultiBasis& e)
{
	m_p0 = e.m_p0;
	m_p1 = e.m_p1;
	m_normal = e.m_normal;
	m_mes = e.m_mes;
	//m_w = e.m_w;
}

void CEdgeMultiBasis::CompNormal()
{
	m_normal.x = m_p0.y - m_p1.y;
	m_normal.y = m_p1.x - m_p0.x;
}

void CEdgeMultiBasis::CompLenght()
{
	m_mes = sqrt((m_p1.x - m_p0.x)*(m_p1.x - m_p0.x) + (m_p1.y - m_p0.y)*(m_p1.y - m_p0.y) + (m_p1.z - m_p0.z)*(m_p1.z - m_p0.z));
}

const double CEdgeMultiBasis::GetShapeFunction(const int k, const Mesh::Point &p) const
{
	// [0;1]
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes };
	const unsigned int n = 10000;
	const double an = 1. / n;
	//vector<double> u(n);
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	if (p.x < m_p0.x || p.x > m_p1.x)
		return 0;
	unsigned int i = 0;
	for (i = 1; i < n ; ++i)
	{
		if (temp.x < (double)i / n)
			break;
	}
	double kek = 0;
	switch (k)
	{
	case 0:
		kek = 1;
		//kek = 1 - ((double)i / n);
		//kek = 1 - ((double)i / n + (double)(i - 1) / n) / 2;
		return kek;
	case 1:
		kek = ((double)i / n);
		//kek = ((double)i / n + (double)(i - 1) / n) / 2;
		return 2 * kek - 1;
	default:
		return -1;
	}
}

const Point CEdgeMultiBasis::GetGradShapeFunction(const int k, const Mesh::Point& p) const
{
	Point temp{ (p.x - m_p0.x) / m_mes,  (p.y - m_p0.y) / m_mes, (p.z - m_p0.z) / m_mes };
	temp.x = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);
	switch (k)
	{
	case 0:
		return Point(0 , 0, 0);
	case 1:
		return Point(1/m_mes, 0, 0);
	}
	return Point(0, 0, 0);
}

const int CEdgeMultiBasis::GetNumberOfShapeFunctions() const
{
	return m_number;
}

const Point CEdgeMultiBasis::GetNormal() const
{
	return m_normal;
}

void CEdgeMultiBasis::ReverseNormal()
{
	m_normal.x = -m_normal.x;
	m_normal.y = -m_normal.y;
	m_normal.z = -m_normal.z;
}

const double CEdgeLinearBasis::GetWeight(const int node, const vector<Point>& pts, const std::function<const double(const Point&)>& func) const
{
    return func(pts[node]);
}


//const CFESolution CEdgeMultiBasis::GetValue(const Point& p) const
//{
//	return m_w[0] * GetShapeFunction(0, p) + m_w[1] * GetShapeFunction(1, p) + m_w[2] * GetShapeFunction(2, p) + m_w[3] * GetShapeFunction(3, p);
//}
//const int CEdgeMultiBasis::SetValue(const int, CSolution* value)
//{
//	return 0;
//}
//const int CEdgeMultiBasis::SetValue(const int n, const CFESolution& value)
//{
//	if (n < m_w.size())
//	{
//		m_w[n] = value;
//		return 0;
//	}
//	return 1;
//}
//const CFESolution CEdgeMultiBasis::GetValue(const int n) const
//{
//	if (n < m_w.size())
//		return m_w[n];
//	return CFESolution(0);
//}
