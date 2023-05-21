#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include "Shape.h"
#include "ShapeFunction.h"
#include <iostream>

namespace corenc
{
	namespace Mesh
	{
		class CNode : public CShape
		{
		public:
			CNode();
			CNode(const CNode&);
			CNode(const int n);
			CNode(const int*n);
			CNode& operator=(const CNode& e)
			{
				m_node = e.m_node;
				return *this;
			}
			friend const bool operator==(const CNode& e1, const CNode& e2)
			{
				if (e1.m_node == e2.m_node)
					return true;
				return false;
			}
			friend std::istream& operator>>(std::istream& is, CNode& e)
			{
				is >> e.m_node;
				--e.m_node;
				return is;
			}
			~CNode() {};
			const int 												GetNode(const int) const;
			const int 												GetNode(const NODES&) const;
			const int												IncreaseOrder() { return 1; };
			const int 												GetNumberOfNodes() const;
			void													SetNode(const int k, const int node);
			const double											Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point												Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const std::vector<double>								Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
		private:
			const int												m_number = 1;
			int														m_node;
		};

		class CNodeBasis : public CShapeFunction<double>
		{
		public:
			CNodeBasis();
			CNodeBasis(const Point*);
			CNodeBasis(const CNodeBasis&e)
			{
				m_p0 = e.m_p0;
				m_normal = e.m_normal;
			};
			CNodeBasis& operator=(const CNodeBasis& e)
			{
				m_p0 = e.m_p0;
				m_normal = e.m_normal;
				return *this;
			}
			~CNodeBasis() {};
			const int												GetNumberOfShapeFunctions() const;
			//const DForm<0>*											GetShapeFunction(const int) const;
			const double											GetShapeFunction(const int, const Point&) const;
			const Point												GetGradShapeFunction(const int, const Point&) const;
			const Point												GetNormal() const;
			void													ReverseNormal();
            const double                                            GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
            {
                return f(verts[node]);
            };
			//const int												SetValue(const int, CSolution* value);
			//const int												SetValue(const int, const CFESolution& value);
			const int												IncreaseOrder() { return 1; };
			//const CFESolution										GetValue(const Point&) const;
			const double											GetMeasure() const { return 0.; };
			//const CFESolution										GetValue(const int) const;
			//const std::function<const DForm<0>*(const Point&)>			GetShapeFunction(const int) const;
		private:
			static const int										m_number = 1;
			Point													m_p0;
			Point													m_normal;
			//CFESolution												m_w;
			//const std::function<const double(const Point&p)>		m_psi[2];
		};
	}

}
#endif /* Node_hpp */
