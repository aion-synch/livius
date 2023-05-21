#include "Mesh1D.h"
#include "../FiniteElements/Node.h"
#include "../FiniteElements/Edge.h"
#include <iostream>
using namespace std;
using namespace corenc;
using namespace Mesh;

CMesh1D::CMesh1D()
{
	m_params.resize(3);
	m_params[0] = 1;
	m_params[1] = 1;
	m_params[2] = 1;
}

CMesh1D::CMesh1D(const std::string& file_name)
{
	m_params.resize(3);
	m_params[0] = 1;
	m_params[1] = 1;
	m_params[2] = 1;
	cout << file_name << endl;
	cout << "Reading the grid from the file..." << endl;
	ifstream in{file_name};
	if(!in.is_open())
	{
		cout << "Error! Couldn't load the grid\n";
		return;
	}
	double a, b, sz;
	unsigned int i, n;
	in >> a >> b >> n;
	sz = (b - a) / n;
	int nodes[2];
	Point pts[2];
	for(i = 0; i < n; ++i)
	{
		nodes[0] = i;
		nodes[1] = i + 1;
		pts[0] = Point(a + i * sz, 0, 0);
		pts[1] = Point(a + (i + 1) * sz, 0, 0);
		m_elems.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, unsigned int, CFESolution>{nodes, pts});
	}
	nodes[0] = 0;
	pts[0] = Point(a,0,0);
	m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{nodes, pts});
	m_bnds[0]->SetNeighbour(0, 0);
	m_bnds[0]->SetNeighbour(1, -1);
	m_points.resize(n + 1);
	m_points[0] = pts[0];
	for(i = 1; i < n; ++i)
	{
		nodes[0] = i;
		pts[0] = Point(a + i * sz, 0, 0);
		m_points[i] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{nodes, pts});
		m_bnds[i]->SetNeighbour(0, i - 1);
		m_bnds[i]->SetNeighbour(1, i);
		/* .-.-.-.-.-.-.-.-.-.
		   0011223344556677889
		 */
	}
	nodes[0] = (int)n;
	pts[0] = Point(b,0,0);
	m_points[n] = pts[0];
	m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{nodes, pts});
	m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
	m_bnds[n]->SetNeighbour(1, -1);
	//cout << "Complete." << endl;
}

CMesh1D::CMesh1D(const std::string& file_name, const std::string& init_file)
{
	m_params.resize(3);
	m_params[0] = 1;
	m_params[1] = 1;
	m_params[2] = 1;
	cout << file_name << endl;
	cout << "Reading the grid from the file..." << endl;
	ifstream in{ file_name };
	if (!in.is_open())
	{
		cout << "Error! Couldn't load the grid\n";
		return;
	}
	double a, b, sz;
	unsigned int i, n;
	in >> a >> b >> n;
	sz = (b - a) / n;
	m_minsize = sz;
	int nodes[2];
	Point pts[2];
	for (i = 0; i < n; ++i)
	{
		nodes[0] = i;
		nodes[1] = i + 1;
		pts[0] = Point(a + i * sz, 0, 0);
		pts[1] = Point(a + (i + 1) * sz, 0, 0);
		m_elems.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, unsigned int, CFESolution>{ nodes, pts });
	}
	nodes[0] = 0;
	pts[0] = Point(a, 0, 0);
	m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
	m_bnds[0]->SetNeighbour(0, 0);
	m_bnds[0]->SetNeighbour(1, -1);
	m_points.resize(n + 1);
	m_points[0] = pts[0];
	for (i = 1; i < n; ++i)
	{
		nodes[0] = i;
		pts[0] = Point(a + i * sz, 0, 0);
		m_points[i] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[i]->SetNeighbour(0, i - 1);
		m_bnds[i]->SetNeighbour(1, i);
		/* .-.-.-.-.-.-.-.-.-.
		0011223344556677889
		*/
	}
	nodes[0] = (int)n;
	pts[0] = Point(b, 0, 0);
	m_points[n] = pts[0];
	m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
	m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
	m_bnds[n]->SetNeighbour(1, -1);
	cout << init_file << endl;
	cout << "Reading the grid from the file..." << endl;
	ifstream init{ init_file };
	if (!init.is_open())
	{
		cout << "Error! Couldn't load the grid\n";
		return;
	}
	init >> n;
	m_solution.resize(n);
	for (auto& it : m_solution)
		init >> it;
	//cout << "Complete." << endl;
}

CMesh1D::CMesh1D(const double a, const double b, const unsigned n, const int order, const std::function<const double(const Point&)>& init_func)
{
	m_params.resize(3);
	m_params[0] = 1;
	m_params[1] = 1;
	m_params[2] = 1;
	double sz;
	unsigned int i;
	sz = (b - a) / n;
	m_minsize = sz;
	int nodes[2];
	CFESolution temp;
	Point pts[2];
	switch (order)
	{
	case 0:
		/*for (i = 0; i < n; ++i)
		{
			nodes[0] = i;
			nodes[1] = i + 1;
			pts[0] = Point(a + i * sz, 0, 0);
			pts[1] = Point(a + (i + 1) * sz, 0, 0);
			//pts[0] = Point(a + i * sz, 0, 0);
			m_elems.push_back(new CFiniteElement<CEdge, CEdgeConstantBasis, unsigned int>{ nodes, pts });
		}

		nodes[0] = 0;
		pts[0] = Point(a, 0, 0);
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int>{ nodes, pts });
		m_bnds[0]->SetNeighbour(0, 0);
		m_bnds[0]->SetNeighbour(1, -1);
		m_points.resize(n);
		m_points[0] = pts[0];
		for (i = 1; i < n; ++i)
		{
			nodes[0] = i;
			pts[0] = Point(a + i * sz, 0, 0);
			m_points[i] = pts[0];
			m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int>{ nodes, pts });
			m_bnds[i]->SetNeighbour(0, i - 1);
			m_bnds[i]->SetNeighbour(1, i);
		}
		nodes[0] = (int)n;
		pts[0] = Point(b, 0, 0);
		m_points[n-1] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int>{ nodes, pts });
		m_bnds[n-1]->SetNeighbour(0, (int)m_elems.size() - 1);
		m_bnds[n-1]->SetNeighbour(1, -1);*/
		for (i = 0; i < n; ++i)
		{
			nodes[0] = i;
			nodes[1] = i + 1;
			pts[0] = Point(a + i * sz, 0, 0);
			pts[1] = Point(a + (i + 1) * sz, 0, 0);
			m_elems.push_back(new CFiniteElement<CEdge, CEdgeConstantBasis, unsigned int, CFESolution>{ nodes, pts });
			//m_elems[i]->SetValue(0, init_func(pts[0]));
			//temp.m_w = init_func(pts[0]);
			//m_elems[i]->SetValue(0, static_cast<CFESolution*>(&temp));
		}
		nodes[0] = 0;
		pts[0] = Point(a, 0, 0);
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[0]->SetNeighbour(0, 0);
		m_bnds[0]->SetNeighbour(1, -1);
		m_points.resize(n + 1);
		m_points[0] = pts[0];
		for (i = 1; i < n; ++i)
		{
			nodes[0] = i;
			pts[0] = Point(a + i * sz, 0, 0);
			m_points[i] = pts[0];
			m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
			m_bnds[i]->SetNeighbour(0, i - 1);
			m_bnds[i]->SetNeighbour(1, i);
			/* .-.-.-.-.-.-.-.-.-.
			0011223344556677889
			*/
		}
		nodes[0] = (int)n;
		pts[0] = Point(b, 0, 0);
		m_points[n] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
		m_bnds[n]->SetNeighbour(1, -1);

		break;
	case 1:
		for (i = 0; i < n; ++i)
		{
			nodes[0] = i;
			nodes[1] = i + 1;
			pts[0] = Point(a + i * sz, 0, 0);
			pts[1] = Point(a + (i + 1) * sz, 0, 0);
			m_elems.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, unsigned int, CFESolution>{ nodes, pts });
		}

		nodes[0] = 0;
		pts[0] = Point(a, 0, 0);
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[0]->SetNeighbour(0, 0);
		m_bnds[0]->SetNeighbour(1, -1);
		m_points.resize(n + 1);
		m_points[0] = pts[0];
		for (i = 1; i < n; ++i)
		{
			nodes[0] = i;
			pts[0] = Point(a + i * sz, 0, 0);
			m_points[i] = pts[0];
			m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
			m_bnds[i]->SetNeighbour(0, i - 1);
			m_bnds[i]->SetNeighbour(1, i);
			/* .-.-.-.-.-.-.-.-.-.
			0011223344556677889
			*/
		}
		nodes[0] = (int)n;
		pts[0] = Point(b, 0, 0);
		m_points[n] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
		m_bnds[n]->SetNeighbour(1, -1);

		break;
	default:
		cout << "WARNING: the order is undefined." << endl;
		break;
	}

	/*for (i = 0; i < n; ++i)
	{
		nodes[0] = i;
		nodes[1] = i + 1;
		pts[0] = Point(a + i * sz, 0, 0);
		pts[1] = Point(a + (i + 1) * sz, 0, 0);
		m_elems.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, unsigned int>{ nodes, pts });
	}*/
	/*nodes[0] = 0;
	pts[0] = Point(a, 0, 0);
	m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int>{ nodes, pts });
	m_bnds[0]->SetNeighbour(0, 0);
	m_bnds[0]->SetNeighbour(1, -1);
	m_points.resize(n + 1);
	m_points[0] = pts[0];
	for (i = 1; i < n; ++i)
	{
		nodes[0] = i;
		pts[0] = Point(a + i * sz, 0, 0);
		m_points[i] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int>{ nodes, pts });
		m_bnds[i]->SetNeighbour(0, i - 1);
		m_bnds[i]->SetNeighbour(1, i);
	}
	nodes[0] = (int)n;
	pts[0] = Point(b, 0, 0);
	m_points[n] = pts[0];
	m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int>{ nodes, pts });
	m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
	m_bnds[n]->SetNeighbour(1, -1);*/
	//CSolution& ke = *m_elems[0]->GetValue(0);
	//CFESolution* ref = dynamic_cast<CFESolution*>(&ke);
	//ref->m_w = 10;
	//m_elems[0]->SetValue(0, ref);
	//CSolution& k2 = *m_elems[0]->GetValue(0);
	//cout << static_cast<CFESolution*>(m_elems[0]->GetValue(0))->m_w << endl;
	//auto& f = static_cast<CFESolution*>(m_elems[0]->GetValue(0));
	//*f = 10.;
	m_solution.resize(m_points.size());
	for (i = 0; i < m_points.size(); ++i)
		m_solution[i] = init_func(m_points[i]);
	//cout << "Complete." << endl;
}

CMesh1D::CMesh1D(const double a, const double b, const unsigned n, const int order, const std::function<const double(const Point&)>& init_func, const std::function<const double(const Point&)>& init_derivative)
{
	m_params.resize(3);
	m_params[0] = 1;
	m_params[1] = 1;
	m_params[2] = 1;
	double sz;
	unsigned int i;
	sz = (b - a) / n;
	m_minsize = sz;
	int nodes[2];
	Point pts[2];
	switch (order)
	{
	case 0:

		for (i = 0; i < n; ++i)
		{
			nodes[0] = i;
			nodes[1] = i + 1;
			pts[0] = Point(a + i * sz, 0, 0);
			pts[1] = Point(a + (i + 1) * sz, 0, 0);
			m_elems.push_back(new CFiniteElement<CEdge, CEdgeHermiteBasis, unsigned int, CFESolution>{ nodes, pts });
			//m_elems[i]->SetValue(0, init_func(pts[0]));
			//m_elems[i]->SetValue(1, init_derivative(pts[0]));
			//m_elems[i]->SetValue(2, init_func(pts[1]));
			//m_elems[i]->SetValue(3, init_derivative(pts[1]));
		}

		nodes[0] = 0;
		pts[0] = Point(a, 0, 0);
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[0]->SetNeighbour(0, 0);
		m_bnds[0]->SetNeighbour(1, -1);
		m_points.resize(n + 1);
		m_points[0] = pts[0];
		for (i = 1; i < n; ++i)
		{
			nodes[0] = i;
			pts[0] = Point(a + i * sz, 0, 0);
			m_points[i] = pts[0];
			m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
			m_bnds[i]->SetNeighbour(0, i - 1);
			m_bnds[i]->SetNeighbour(1, i);
		}
		nodes[0] = (int)n;
		pts[0] = Point(b, 0, 0);
		m_points[n] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
		m_bnds[n]->SetNeighbour(1, -1);

		break;
	case 1:
		for (i = 0; i < n; ++i)
		{
			nodes[0] = i;
			nodes[1] = i + 1;
			pts[0] = Point(a + i * sz, 0, 0);
			pts[1] = Point(a + (i + 1) * sz, 0, 0);
			m_elems.push_back(new CFiniteElement<CEdge, CEdgeLinearBasis, unsigned int, CFESolution>{ nodes, pts });
		}

		nodes[0] = 0;
		pts[0] = Point(a, 0, 0);
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[0]->SetNeighbour(0, 0);
		m_bnds[0]->SetNeighbour(1, -1);
		m_points.resize(n + 1);
		m_points[0] = pts[0];
		for (i = 1; i < n; ++i)
		{
			nodes[0] = i;
			pts[0] = Point(a + i * sz, 0, 0);
			m_points[i] = pts[0];
			m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
			m_bnds[i]->SetNeighbour(0, i - 1);
			m_bnds[i]->SetNeighbour(1, i);
			/* .-.-.-.-.-.-.-.-.-.
			0011223344556677889
			*/
		}
		nodes[0] = (int)n;
		pts[0] = Point(b, 0, 0);
		m_points[n] = pts[0];
		m_bnds.push_back(new CFiniteElement<CNode, CNodeBasis, unsigned int, CFESolution>{ nodes, pts });
		m_bnds[n]->SetNeighbour(0, (int)m_elems.size() - 1);
		m_bnds[n]->SetNeighbour(1, -1);

		break;
	default:
		cout << "WARNING: the order is undefined." << endl;
		break;
	}

	m_solution.resize(2 * m_points.size());
	for (i = 0; i < m_points.size(); ++i)
	{
		m_solution[2 * i] = init_func(m_points[i]);
	}
	//cout << "Complete." << endl;
}

CMesh1D::CMesh1D(const CMesh1D& m)
{
	auto sz{m.m_elems.size()};
	auto sz2{m.m_bnds.size()};
	auto sz3{m.m_points.size()};
	unsigned int i{0};
	m_elems.resize(sz);
	m_bnds.resize(sz2);
	m_points.resize(sz3);
	for(; i < sz; ++i)
		m_elems[i] = m.m_elems[i]->Clone();
	for(i = 0; i < sz2; ++i)
		m_bnds[i] = m.m_bnds[i]->Clone();
	for(i = 0; i < sz3; ++i)
		m_points[i] = m.m_points[i];
}

CMesh1D::~CMesh1D()
{
	if(m_elems.size())
	{
		auto sz{m_elems.size()};
		for(auto i{0}; i < sz; ++i)
			delete m_elems[i];
		vector<CElement<CFESolution>*>().swap(m_elems);
	}
	if(m_bnds.size())
	{
		auto sz{m_bnds.size()};
		for(auto i{0}; i < sz; ++i)
			delete m_bnds[i];
		vector<CElement<CFESolution>*>().swap(m_bnds);
	}
}

const unsigned int CMesh1D::GetNumberOfElements() const
{
	return (unsigned int)m_elems.size();
}

const unsigned int CMesh1D::GetNumberOfNodes() const
{
	return (unsigned int)m_points.size();
}

const unsigned int CMesh1D::GetNumberOfBoundaries() const
{
	return (unsigned int)m_bnds.size();
}

const int CMesh1D::FindElement(const Point& test) const
{
	auto sz{m_elems.size()};
	Point p[4];
	const double temp{ sqrt(test.x *test.x + test.y * test.y) };
	for (unsigned int i = 0; i < sz; ++i)
	{
		p[0].x = sqrt(m_points[m_elems[i]->GetNode(0)].x * m_points[m_elems[i]->GetNode(0)].x + m_points[m_elems[i]->GetNode(0)].y*m_points[m_elems[i]->GetNode(0)].y);
		p[1].x = sqrt(m_points[m_elems[i]->GetNode(1)].x * m_points[m_elems[i]->GetNode(1)].x + m_points[m_elems[i]->GetNode(1)].y*m_points[m_elems[i]->GetNode(1)].y);
		if (p[0].x < p[1].x)
		{
			if (temp >= p[0].x && temp <= p[1].x)
				return i;
		}
		else
			if (temp <= p[0].x && temp >= p[1].x)
				return i;
		if (fabs(temp - p[0].x) < 1e-13)
			return i;
		if (fabs(temp - p[1].x) < 1e-13)
			return i;
	}
	return -1;
}

const CElement<CFESolution>* CMesh1D::GetElement(const unsigned int n) const
{
	if(n < m_elems.size())
		return m_elems[n];
	cout << "Wrong number of the element." << endl;
	return nullptr;
}

const CElement<CFESolution>* CMesh1D::GetBoundary(const unsigned int n) const
{
	if(n < m_bnds.size())
		return m_bnds[n];
	cout << "Wrong number of the boundary." << endl;
	return nullptr;
}

const Point CMesh1D::GetNode(const unsigned int n) const
{
	if(n < m_points.size())
		return m_points[n];
	cout << "Wrong number of the node." << endl;
	return Point(0,0,0);
}

const int CMesh1D::updateSolution(const vector<double>& new_solution)
{
	if (new_solution.size() != m_solution.size())
		return 1;
	m_solution = new_solution;
	return 0;
}

const int CMesh1D::updateSolution(const unsigned int element, const unsigned int node, CSolution* value)
{
	//m_elems[element]->SetValue(node, value);
	return 0;
}

const double CMesh1D::getSolution(const unsigned int element, const unsigned int node) const
{
	//return m_elems[element]->GetValue(node);
	return m_solution[m_elems[element]->GetNode(node)];
}

const int CMesh1D::updateSolution(const unsigned int element, const unsigned int node, const double value)
{
	m_solution[element + node] = value;
	//m_elems[element]->SetValue(node, value);
	//m_solution[element + node] = 1./(double)m_elems[element]->GetValue(node);
	return 0;
}

const int CMesh1D::updateSolution(const unsigned int node, const double value)
{
	if (node < m_solution.size())
	{
		m_solution[node] = value;
		return 0;
	}
	cout << "WARNING: UPDATE WRONG NODE! CMESH1D::updateSolution(node, value)" << endl;
	return 1;
}

const double CMesh1D::getParameter(Parameters param, const unsigned int, const Point&) const
{
	switch (param)
	{
	case corenc::Parameters::DIFFUSION:
		return m_params[0];
		break;
	case corenc::Parameters::MASS:
		return m_params[1];
		break;
	case corenc::Parameters::ADVECTION:
		return m_params[2];
		break;
	default:
		break;
	}
	return 0;
}

const double CMesh1D::getParameter(Parameters param, const unsigned int, const int) const
{
	switch (param)
	{
	case corenc::Parameters::DIFFUSION:
		return m_params[0];
		break;
	case corenc::Parameters::MASS:
		return m_params[1];
		break;
	case corenc::Parameters::ADVECTION:
		return m_params[2];
		break;
	default:
		break;
	}
	return 0;
}

const int CMesh1D::setParameter(Parameters param, const double value, const unsigned int)
{
	switch (param)
	{
	case corenc::Parameters::DIFFUSION:
		m_params[0] = value;
		break;
	case corenc::Parameters::MASS:
		m_params[1] = value;
		break;
	case corenc::Parameters::ADVECTION:
		m_params[2] = value;
		break;
	default:
		break;
	}
	return 0;
}
