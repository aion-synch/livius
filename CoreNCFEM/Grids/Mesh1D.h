#ifndef CORENC_Mesh1D_hpp
#define CORENC_Mesh1D_hpp

#include <stdio.h>
#include "../Mesh.h"

#include "../Point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <typeinfo>
namespace corenc
{
	namespace Mesh
	{
		
		class CMesh1D : public CMesh<CFESolution>
		{
		public:
			CMesh1D();
			CMesh1D(const std::string& domain_name);
			CMesh1D(const std::string& domain_file, const std::string& init_file);
			CMesh1D(const double x0, const double x1, const unsigned n, const int order, const std::function<const double(const Point&)>& init_func);
			CMesh1D(const double x0, const double x1, const unsigned n, const int order, const std::function<const double(const Point&)>& init_func, const std::function<const double(const Point&)>& init_derivative);
			CMesh1D(const CMesh1D&);
			CMesh1D& operator=(const CMesh1D& m)
			{
				auto sz = m.m_elems.size();
				m_elems.resize(sz);
				for (int i = 0; i < sz; ++i)
					m_elems[i] = m.m_elems[i]->Clone();
			}
			const unsigned int							GetNumberOfElements() const;
			const unsigned int							GetNumberOfNodes() const;
			const unsigned int							GetNumberOfBoundaries() const;
			const int									FindElement(const Point&) const;
			const Point									GetNode(const unsigned int) const;
			const CElement<CFESolution>*				GetElement(const unsigned int) const;
			const CElement<CFESolution>*				GetBoundary(const unsigned int) const;
			const double								getSolution(const unsigned int element, const unsigned int node) const;
			const double								getParameter(Parameters, const unsigned int, const Point& p) const;
			const double								getParameter(Parameters, const unsigned int, const int) const;
			const std::vector<double>					getSolution() const { return m_solution; };
			const int									updateSolution(const std::vector<double>& new_solution);
			const int									updateSolution(const unsigned int element, const unsigned int node, const double value);
			const int									updateSolution(const unsigned int element, const unsigned int node, CSolution* value);
			const int									updateSolution(const unsigned int node, const double value);
			const int									setParameter(Parameters, const double, const unsigned int);
			const double								getMinSize() const { return m_minsize; };
			~CMesh1D();
		private:
			std::vector<CElement<CFESolution>*>			m_elems;
			std::vector<CElement<CFESolution>*>			m_bnds;
			std::vector<Point>							m_points;
			std::vector<double>							m_solution;
			std::vector<int>							m_nums;
			std::vector<double>							m_params;
			double										m_minsize{0.};
		public:
			auto										GetElements() -> decltype(m_elems) { return m_elems; };
			auto										GetBoundary() -> decltype(m_bnds) { return m_bnds; };
		};
	}
}
#endif /* CORENC_Mesh1D_hpp */
