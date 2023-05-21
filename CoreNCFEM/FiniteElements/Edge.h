#ifndef Edge_hpp
#define Edge_hpp

#include <stdio.h>
#include "Shape.h"
#include "ShapeFunction.h"
#include <iostream>
#include "../FESolution.h"
namespace corenc
{
	namespace Mesh
	{
		class CEdge : public CShape
		{
		public:
			CEdge();
			CEdge(const CEdge&);
			CEdge(const int n1, const int n2);
			CEdge(const int*);
			CEdge& operator=(const CEdge& e)
			{
				m_nodes = e.m_nodes;
				m_number = e.m_number;
				return *this;
			}
			friend const bool operator==(const CEdge& e1, const CEdge& e2)
			{
				if (e1.m_nodes[0] == e2.m_nodes[0])
					if (e1.m_nodes[1] == e2.m_nodes[1])
						return true;
				if (e1.m_nodes[1] == e2.m_nodes[0])
					if (e1.m_nodes[0] == e2.m_nodes[1])
						return true;
				return false;
			}
			friend std::istream& operator>>(std::istream& is, CEdge& e)
			{
				is >> e.m_nodes[0] >> e.m_nodes[1];
				--e.m_nodes[0]; --e.m_nodes[1];
				return is;
			}
			~CEdge() {};
			const int 												GetNode(const int) const;
			const int 												GetNode(const NODES&) const;
			const int 												GetNumberOfNodes() const;
			void													SetNode(const int k, const int node);
			const double											Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point												Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const int												IncreaseOrder();
			const std::vector<double>								Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
		private:
			int														m_number;
			std::vector<int>										m_nodes;
		};

		class CEdgeLinearBasis : public CShapeFunction<double>
		{
		public:
			CEdgeLinearBasis();
			CEdgeLinearBasis(const Point&, const Point&);
			CEdgeLinearBasis(const Point*);
			CEdgeLinearBasis(const CEdgeLinearBasis&);
			CEdgeLinearBasis& operator=(const CEdgeLinearBasis& e)
			{
				m_number = e.m_number;
				m_p0 = e.m_p0;
				m_p1 = e.m_p1;
				m_normal = e.m_normal;
				m_mes = e.m_mes;
				return *this;
			}
			~CEdgeLinearBasis() {};
			const int												GetNumberOfShapeFunctions() const;
			//const DForm<0>*											GetShapeFunction(const int) const;
			const double											GetShapeFunction(const int, const Point&) const;
			const Point												GetGradShapeFunction(const int, const Point&) const;
			const Point												GetNormal() const;
			void													ReverseNormal();
			//const int												SetValue(const int, CSolution* value);
			//const int												SetValue(const int, const CFESolution& value);
			//const CFESolution										GetValue(const Point&) const;
			//const CFESolution										GetValue(const int) const;
			const int												IncreaseOrder();
			const double											GetMeasure() const { return m_mes; };
            const double                                            GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
			//CSolution*										GetValue(const unsigned int);
			//const std::function<const DForm<0>*(const Point&)>			GetShapeFunction(const int) const;
		private:
			int														m_number;
			Point													m_p0, m_p1;
			Point													m_normal;
			double													m_mes;
			void													CompNormal();
			void													CompLenght();
			//std::vector<double>										m_w[2];
			//std::vector<CFESolution>								m_w;

			//const std::function<const double(const Point&p)>		m_psi[2];
		};



		class CEdgeConstantBasis : public CShapeFunction<double>
		{
		public:
			CEdgeConstantBasis();
			CEdgeConstantBasis(const Point&, const Point&);
			CEdgeConstantBasis(const Point*);
			CEdgeConstantBasis(const CEdgeConstantBasis&);
			CEdgeConstantBasis& operator=(const CEdgeConstantBasis& e)
			{
				m_p0 = e.m_p0;
				m_p1 = e.m_p1;
				m_normal = e.m_normal;
				m_mes = e.m_mes;
				return *this;
			}
			~CEdgeConstantBasis() {};
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
			//const double											GetValue(const Point&) const;
			//const int												SetValue(const unsigned int, const double& value);
			//const double											GetValue(const unsigned int) const;
			//const int												SetValue(const int, CSolution* value);
			//const int												SetValue(const int, const CFESolution& value);
			const int												IncreaseOrder() { return 1; };
			//const CFESolution										GetValue(const Point&) const;
			//CSolution*										GetValue(const unsigned int);
			//const CFESolution										GetValue(const int) const;
			const double											GetMeasure() const { return 0.; };
			//const std::function<const DForm<0>*(const Point&)>			GetShapeFunction(const int) const;
		private:
			static const int										m_number = 1;
			Point													m_p0;
			Point													m_p1;
			Point													m_normal;
			double													m_mes;
			void													CompNormal();
			void													CompLenght();
			//double													m_w;
			//CFESolution												m_w;
			//const std::function<const double(const Point&p)>		m_psi[2];
		};

		class CEdgeMultiBasis : public CShapeFunction<double>
		{
		public:
			CEdgeMultiBasis();
			CEdgeMultiBasis(const Point&, const Point&);
			CEdgeMultiBasis(const Point*);
			CEdgeMultiBasis(const CEdgeMultiBasis&);
			CEdgeMultiBasis& operator=(const CEdgeMultiBasis& e)
			{
				m_p0 = e.m_p0;
				m_p1 = e.m_p1;
				m_normal = e.m_normal;
				m_mes = e.m_mes;
				//m_w = e.m_w;
				return *this;
			}
			~CEdgeMultiBasis() {};
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
			//const int												SetValue(const  int, CSolution* value);
			const int												IncreaseOrder() { return 1; };
			//const int												SetValue(const int, const CFESolution& value);
			//const CFESolution										GetValue(const Point&) const;
			//const CFESolution										GetValue(const int) const;
            const double											GetMeasure() const { return 0.; };
			//const std::function<const DForm<0>*(const Point&)>			GetShapeFunction(const int) const;
		private:
			static const int										m_number = 2;
			Point													m_p0;
			Point													m_p1;
			Point													m_normal;
			double													m_mes;
			void													CompNormal();
			void													CompLenght();
			//std::vector<double>										m_w[m_number];
			//std::vector<CFESolution>								m_w;
			//const std::function<const double(const Point&p)>		m_psi[2];
		};

		class CEdgeHermiteBasis : public CShapeFunction<double>
		{
		public:
			CEdgeHermiteBasis();
			CEdgeHermiteBasis(const Point&, const Point&);
			CEdgeHermiteBasis(const Point*);
			CEdgeHermiteBasis(const CEdgeHermiteBasis&);
			CEdgeHermiteBasis& operator=(const CEdgeHermiteBasis& e)
			{
				m_p0 = e.m_p0;
				m_p1 = e.m_p1;
				m_normal = e.m_normal;
				m_mes = e.m_mes;
				//m_w = e.m_w;
				return *this;
			}
			~CEdgeHermiteBasis() {};
			const int												GetNumberOfShapeFunctions() const;
			//const DForm<0>*											GetShapeFunction(const int) const;
			const double											GetShapeFunction(const int, const Point&) const;
			const Point												GetGradShapeFunction(const int, const Point&) const;
			const int												IncreaseOrder() { return 1; };
			const Point												GetNormal() const;
			void													ReverseNormal();
            const double                                            GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
            {
                return f(verts[node]);
            };
			//const int												SetValue(const int, CSolution* value);
			//const int												SetValue(const int, const CFESolution& value);
			//const CFESolution										GetValue(const Point&) const;
			//const CFESolution										GetValue(const int) const;
			const double											GetMeasure() const { return 0.; };
			//const std::function<const DForm<0>*(const Point&)>			GetShapeFunction(const int) const;
		private:
			static const int										m_number = 4;
			Point													m_p0;
			Point													m_p1;
			Point													m_normal;
			double													m_mes;
			void													CompNormal();
			void													CompLenght();
			//std::vector<double>										m_w[m_number];
			//std::vector<CFESolution>								m_w;
			//const std::function<const double(const Point&p)>		m_psi[2];
		};

		class CEdge2ndBasis : public CShapeFunction<double>
		{
		public:
			CEdge2ndBasis();
			CEdge2ndBasis(const Point&, const Point&);
			CEdge2ndBasis(const Point*);
			CEdge2ndBasis(const CEdge2ndBasis&);
			CEdge2ndBasis& operator=(const CEdge2ndBasis& e)
			{
				m_p0 = e.m_p0;
				m_p1 = e.m_p1;
				m_normal = e.m_normal;
				m_mes = e.m_mes;
				return *this;
			}
			~CEdge2ndBasis() {};
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
			const int												IncreaseOrder() { return 1; };
			//const int												SetValue(const int, CSolution* value);
			//const int												SetValue(const int, const CFESolution& value);
			//const CFESolution										GetValue(const Point&) const;
			const double											GetMeasure() const { return 0.; };
			//const CFESolution										GetValue(const int) const;
			//const std::function<const DForm<0>*(const Point&)>			GetShapeFunction(const int) const;
		private:
			static const int										m_number = 3;
			Point													m_p0, m_p1;
			Point													m_normal;
			double													m_mes;
			void													CompNormal();
			void													CompLenght();
			//std::vector<CFESolution>								m_w;
			//const std::function<const double(const Point&p)>		m_psi[2];
		};
	}
}

#endif /* Edge_hpp */
