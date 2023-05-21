#ifndef CORENC_MESH_TRIANGLELINEAR_H_
#define CORENC_MESH_TRIANGLELINEAR_H_

#include <stdio.h>
#include "Shape.h"
#include "ShapeFunction.h"
#include <iostream>
namespace corenc
{
	namespace Mesh
	{
		class CTriangleLinear : public CShape
		{
		public:
			CTriangleLinear();
			CTriangleLinear(const int n1, const int n2, const int n3);
			CTriangleLinear(const int n1, const int n2, const int n3, const int e1, const int e2, const int e3);
			CTriangleLinear(const int*);
			CTriangleLinear(const int*, const int*);
			CTriangleLinear(const CTriangleLinear&);
			CTriangleLinear& operator=(const CTriangleLinear& t)
			{
				m_nodes[0] = t.m_nodes[0];
				m_nodes[1] = t.m_nodes[1];
				m_nodes[2] = t.m_nodes[2];
				m_edges[0] = t.m_edges[0];
				m_edges[1] = t.m_edges[1];
				m_edges[2] = t.m_edges[2];
				return *this;
			}
			const bool	  operator==(const CTriangleLinear& t)
			{
				for (unsigned int i = 0; i < 3; ++i)
					if (m_nodes[i] == t.m_nodes[0])
						for (unsigned int j = 0; j < 3; ++j)
							if (m_nodes[j] == t.m_nodes[1])
								for (unsigned int k = 0; k < 3; ++k)
									if (m_nodes[k] == t.m_nodes[2])
										return true;
				return false;
			}
			std::istream& operator>>(std::istream& is)
			{
				is >> m_nodes[0] >> m_nodes[1] >> m_nodes[2];
				return is;
			}
			~CTriangleLinear() {};
			const int                                         	GetNode(const int) const;
			const int                                         	GetNode(const NODES&) const;
			const int                                           GetEdge(const int) const;
			const int                                           GetFacet(const int) const;
			const int			                                GetNumberOfNodes() const;
			const int			                                GetNumberOfEdges() const;
			const int											GetNumberOfFacets() const;
			const double										Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point											Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const std::vector<double>							Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
			void                                                SetNode(const int k, const int node);
			const int											IncreaseOrder() { return 1; };
			void                                                SetEdge(const int k, const int edge);
			void                                                SetFacet(const int k, const int facet);
		private:
			int                                                 m_nodes[3];
			int                                                 m_edges[3];
		};

		class CTriangleLinearBasis : public CShapeFunction<double>
		{
		public:
			CTriangleLinearBasis();
			CTriangleLinearBasis(const Point&, const Point&, const Point&);
			CTriangleLinearBasis(const Point*);
			CTriangleLinearBasis(const CTriangleLinearBasis&);
			CTriangleLinearBasis& operator=(const CTriangleLinearBasis& t)
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
				return *this;
			}
			~CTriangleLinearBasis() {};
			const int				                            GetNumberOfShapeFunctions() const;
			//const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
			const double				                        GetShapeFunction(const int, const Point&) const;
			const Point											GetGradShapeFunction(const int, const Point&) const;
			const Point											GetNormal() const;
			void												ReverseNormal();
			const double										GetValue(const Point&) const;
			const int											IncreaseOrder() { return 1; };
			//const int											SetValue(const int, CSolution* value);
			//CSolution*										GetValue(const unsigned int);
			//const CFESolution									GetValue(const int) const;
			const double										GetMeasure() const { return m_det; };
			//const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
		private:
			static const int                                    m_number = 3;
			double                                              m_alpha[3][3];
			Point                                               m_normal;
			double                                              m_det;
			const double                                        L(const int, const Point&) const;
			void                                                compD(const Point&, const Point&, const Point&);
			void                                                compAlpha(const Point&, const Point&, const Point&);
			void                                                compNormal(const Point&, const Point&, const Point&);
			//std::vector<double>									m_w[m_number];
			//std::vector<CFESolution>							m_w{3};
		};

		class CTriangleBasis : public CShapeFunction<double>
		{
		public:
			CTriangleBasis();
			CTriangleBasis(const Point&, const Point&, const Point&, const int order);
			CTriangleBasis(const Point*, const int order);
			CTriangleBasis(const CTriangleBasis&);
			CTriangleBasis& operator=(const CTriangleBasis& t)
			{
				m_normal = t.m_normal;
				m_det = t.m_det;
				m_order = t.m_order;
				m_number = t.m_number;
				m_alpha[0][0] = t.m_alpha[0][0];
				m_alpha[0][1] = t.m_alpha[0][1];
				m_alpha[0][2] = t.m_alpha[0][2];

				m_alpha[1][0] = t.m_alpha[1][0];
				m_alpha[1][1] = t.m_alpha[1][1];
				m_alpha[1][2] = t.m_alpha[1][2];

				m_alpha[2][0] = t.m_alpha[2][0];
				m_alpha[2][1] = t.m_alpha[2][1];
				m_alpha[2][2] = t.m_alpha[2][2];
				return *this;
			}
			~CTriangleBasis() {};
			const int				                            GetNumberOfShapeFunctions() const;
			//const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
			const double				                        GetShapeFunction(const int, const Point&) const;
			const Point											GetGradShapeFunction(const int, const Point&) const;
			const Point											GetNormal() const;
			void												ReverseNormal();
			const double										GetValue(const Point&) const;
			//const int											SetValue(const int, CSolution* value);
			//CSolution*										GetValue(const unsigned int);
			//const CFESolution									GetValue(const int) const;
			//const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
		private:
			int													m_number;
			int													m_order;
			double                                              m_alpha[3][3];
			Point                                               m_normal;
			double                                              m_det;
			const double                                        L(const int, const Point&) const;
			void                                                compD(const Point&, const Point&, const Point&);
			void                                                compAlpha(const Point&, const Point&, const Point&);
			void                                                compNormal(const Point&, const Point&, const Point&);
			//std::vector<double>									m_w[m_number];
			//std::vector<CFESolution>							m_w{ 3 };
		};
	}
}
#endif // CORENC_MESH_TRIANGLELINEAR_H_
