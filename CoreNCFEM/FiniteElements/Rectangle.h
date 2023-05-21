#ifndef CORENC_MESH_RECTANGLE_H_
#define CORENC_MESH_RECTANGLE_H_

#include <stdio.h>
#include "Shape.h"
#include "ShapeFunction.h"
#include <iostream>
namespace corenc
{
	namespace Mesh
	{
		class CRectangle : public CShape
		{
		public:
			CRectangle();
			CRectangle(const int n1, const int n2, const int n3, const int n4, const int order);
			CRectangle(const int n1, const int n2, const int n3, const int n4, const int e1, const int e2, const int e3, const int e4, const int order);
			CRectangle(const int*, const int order);
			CRectangle(const int*, const int*, const int order);
			CRectangle(const CRectangle&);
			CRectangle& operator=(const CRectangle& t)
			{
				m_nodes = t.m_nodes;
				m_edges[0] = t.m_edges[0];
				m_edges[1] = t.m_edges[1];
				m_edges[2] = t.m_edges[2];
				m_edges[3] = t.m_edges[3];
				m_number = t.m_number;
				m_order = t.m_order;
                m_px = t.m_px;
                m_py = t.m_py;
				return *this;
			}
			const bool	  operator==(const CRectangle& t)
			{
				for (unsigned int i = 0; i < 4; ++i)
					if (m_nodes[i] == t.m_nodes[0])
						for (unsigned int j = 0; j < 4; ++j)
							if (m_nodes[j] == t.m_nodes[1])
								for (unsigned int k = 0; k < 4; ++k)
									if (m_nodes[k] == t.m_nodes[2])
										for (unsigned int l = 0; l < 4; ++l)
											if (m_nodes[l] == t.m_nodes[3])
										return true;
				return false;
			}
			std::istream& operator>>(std::istream& is)
			{
				is >> m_nodes[0] >> m_nodes[1] >> m_nodes[2] >> m_nodes[3];
				return is;
			}
			~CRectangle() {};
			const int                                         	GetNode(const int) const;
			const int                                         	GetNode(const NODES&) const;
			const int                                           GetEdge(const int) const;
			const int                                           GetFacet(const int) const;
			const int			                                GetNumberOfNodes() const;
			const int				                            GetNumberOfEdges() const;
			const int											GetNumberOfFacets() const;
			const double										Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point											Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const std::vector<double>							Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
			void                                                SetNode(const int k, const int node);
			const int											IncreaseOrder();
            const int                                           SetOrder(const int px, const int py);
			void                                                SetEdge(const int k, const int edge);
			void                                                SetFacet(const int k, const int facet);
		private:
			std::vector<int>                                    m_nodes;
			int				                                    m_edges[4];
			int													m_order;
			int													m_number;
            int                                                 m_px, m_py;
		};

		class CRectangleBasis : public CShapeFunction<double>
		{
		public:
			CRectangleBasis();
			CRectangleBasis(const Point&, const Point&, const Point&, const Point&, const int order);
			CRectangleBasis(const Point*, const int order);
			CRectangleBasis(const CRectangleBasis&);
			CRectangleBasis& operator=(const CRectangleBasis& t)
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
				return *this;
			}
			~CRectangleBasis() {};
			const int				                            GetNumberOfShapeFunctions() const;
			//const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
			const double				                        GetShapeFunction(const int, const Point&) const;
			const Point											GetGradShapeFunction(const int, const Point&) const;
			const Point											GetNormal() const;
			void												ReverseNormal();
			const double										GetValue(const Point&) const;
			const int											IncreaseOrder();
			//const int											SetValue(const int, CSolution* value);
			//CSolution*										GetValue(const unsigned int);
			//const CFESolution									GetValue(const int) const;
			const double										GetMeasure() const { return m_det; };
			const double										GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
			//const unsigned int									GetOrder() const;
			//const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
		private:
			int													m_number;
			int													m_order;
			int													m_1dorder;
			Point                                               m_normal;
			std::vector<Mesh::Point>							m_points;
			double                                              m_det;
			double												m_hx, m_hy;
			const double										m_x1(const double) const;
			const double										m_x2(const double) const;
			const double										m_y1(const double) const;
			const double										m_y2(const double) const;
			void                                                compD(const Point&, const Point&, const Point&);
			void                                                compNormal(const Point&, const Point&, const Point&);
			const int											createS();
			//std::vector<double>									m_w[m_number];
			//std::vector<CFESolution>							m_w{ 4 };
			int													m_s;
			int													m_sp;
		};


        class CRectangleHBasis : public CShapeFunction<double>
        {
        public:
            CRectangleHBasis();
            CRectangleHBasis(const Point&, const Point&, const Point&, const Point&, const int order);
            CRectangleHBasis(const Point&, const Point&, const Point&, const Point&, const int px, const int py);
            CRectangleHBasis(const Point*, const int order);
            CRectangleHBasis(const Point*, const int px, const int py);
            CRectangleHBasis(const CRectangleHBasis&);
            CRectangleHBasis& operator=(const CRectangleHBasis& t)
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
                m_px = t.m_px;
                m_py = t.m_py;
                return *this;
            }
            ~CRectangleHBasis() {};
            const int				                            GetNumberOfShapeFunctions() const;
            //const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
            const double				                        GetShapeFunction(const int, const Point&) const;
            const Point											GetGradShapeFunction(const int, const Point&) const;
            const Point											GetNormal() const;
            void												ReverseNormal();
            const double										GetValue(const Point&) const;
            const int											IncreaseOrder();
            const int                                           SetOrder(const int px, const int py);
            //const int											SetValue(const int, CSolution* value);
            //CSolution*										GetValue(const unsigned int);
            //const CFESolution									GetValue(const int) const;
            const double										GetMeasure() const { return m_det; };
            const double										GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
            //const unsigned int									GetOrder() const;
            //const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
        private:
            int													m_number;
            int													m_order;
            int													m_1dorder;
            Point                                               m_normal;
            std::vector<Mesh::Point>							m_points;
            double                                              m_det;
            double												m_hx, m_hy;
            int                                                 m_px, m_py;
            const double										m_x1(const double) const;
            const double										m_x2(const double) const;
            const double										m_y1(const double) const;
            const double										m_y2(const double) const;
            const double                                        m_xi(const double, const int n) const;
            const double                                        m_dxi(const double, const int n) const;
            void                                                compD(const Point&, const Point&, const Point&);
            void                                                compNormal(const Point&, const Point&, const Point&);
            const int											createS();
            //std::vector<double>									m_w[m_number];
            //std::vector<CFESolution>							m_w{ 4 };
            int													m_s;
            int													m_sp;
        };

        class CRectangleBasis2x : public CShapeFunction<double>
        {
        public:
            CRectangleBasis2x();
            CRectangleBasis2x(const Point&, const Point&, const Point&, const Point&, const int order);
            CRectangleBasis2x(const Point*, const int order);
            CRectangleBasis2x(const CRectangleBasis2x&);
            CRectangleBasis2x& operator=(const CRectangleBasis2x& t)
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
                return *this;
            }
            ~CRectangleBasis2x() {};
            const int				                            GetNumberOfShapeFunctions() const;
            //const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
            const double				                        GetShapeFunction(const int, const Point&) const;
            const Point											GetGradShapeFunction(const int, const Point&) const;
            const Point											GetNormal() const;
            void												ReverseNormal();
            const double										GetValue(const Point&) const;
            const int											IncreaseOrder();
            //const int											SetValue(const int, CSolution* value);
            //CSolution*										GetValue(const unsigned int);
            //const CFESolution									GetValue(const int) const;
            const double										GetMeasure() const { return m_det; };
            const double										GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
            //const unsigned int									GetOrder() const;
            //const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
        private:
            int													m_number;
            int													m_order;
            int													m_1dorder;
            Point                                               m_normal;
            std::vector<Mesh::Point>							m_points;
            double                                              m_det;
            double												m_hx, m_hy;
            const double										m_x1(const double) const;
            const double										m_x2(const double) const;
            const double										m_x3(const double) const;
            const double                                        m_dx3(const double) const;
            const double										m_y1(const double) const;
            const double										m_y2(const double) const;
            void                                                compD(const Point&, const Point&, const Point&);
            void                                                compNormal(const Point&, const Point&, const Point&);
            const int											createS();
            //std::vector<double>									m_w[m_number];
            //std::vector<CFESolution>							m_w{ 4 };
            int													m_s;
            int													m_sp;
        };

        class CRectangleBasis2y : public CShapeFunction<double>
        {
        public:
            CRectangleBasis2y();
            CRectangleBasis2y(const Point&, const Point&, const Point&, const Point&, const int order);
            CRectangleBasis2y(const Point*, const int order);
            CRectangleBasis2y(const CRectangleBasis2y&);
            CRectangleBasis2y& operator=(const CRectangleBasis2y& t)
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
                return *this;
            }
            ~CRectangleBasis2y() {};
            const int				                            GetNumberOfShapeFunctions() const;
            //const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
            const double				                        GetShapeFunction(const int, const Point&) const;
            const Point											GetGradShapeFunction(const int, const Point&) const;
            const Point											GetNormal() const;
            void												ReverseNormal();
            const double										GetValue(const Point&) const;
            const int											IncreaseOrder();
            //const int											SetValue(const int, CSolution* value);
            //CSolution*										GetValue(const unsigned int);
            //const CFESolution									GetValue(const int) const;
            const double										GetMeasure() const { return m_det; };
            const double										GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
            //const unsigned int									GetOrder() const;
            //const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
        private:
            int													m_number;
            int													m_order;
            int													m_1dorder;
            Point                                               m_normal;
            std::vector<Mesh::Point>							m_points;
            double                                              m_det;
            double												m_hx, m_hy;
            const double										m_x1(const double) const;
            const double										m_x2(const double) const;
            const double										m_y3(const double) const;
            const double                                        m_dy3(const double) const;
            const double										m_y1(const double) const;
            const double										m_y2(const double) const;
            void                                                compD(const Point&, const Point&, const Point&);
            void                                                compNormal(const Point&, const Point&, const Point&);
            const int											createS();
            //std::vector<double>									m_w[m_number];
            //std::vector<CFESolution>							m_w{ 4 };
            int													m_s;
            int													m_sp;
        };

        class CRectangleBasis2 : public CShapeFunction<double>
        {
        public:
            CRectangleBasis2();
            CRectangleBasis2(const Point&, const Point&, const Point&, const Point&, const int order);
            CRectangleBasis2(const Point*, const int order);
            CRectangleBasis2(const CRectangleBasis2&);
            CRectangleBasis2& operator=(const CRectangleBasis2& t)
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
                return *this;
            }
            ~CRectangleBasis2() {};
            const int				                            GetNumberOfShapeFunctions() const;
            //const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
            const double				                        GetShapeFunction(const int, const Point&) const;
            const Point											GetGradShapeFunction(const int, const Point&) const;
            const Point											GetNormal() const;
            void												ReverseNormal();
            const double										GetValue(const Point&) const;
            const int											IncreaseOrder();
            //const int											SetValue(const int, CSolution* value);
            //CSolution*										GetValue(const unsigned int);
            //const CFESolution									GetValue(const int) const;
            const double										GetMeasure() const { return m_det; };
            const double										GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
            //const unsigned int									GetOrder() const;
            //const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
        private:
            int													m_number;
            int													m_order;
            int													m_1dorder;
            Point                                               m_normal;
            std::vector<Mesh::Point>							m_points;
            double                                              m_det;
            double												m_hx, m_hy;
            const double										m_x1(const double) const;
            const double										m_x2(const double) const;
            const double										m_x3(const double) const;
            const double                                        m_dx3(const double) const;
            const double										m_y1(const double) const;
            const double										m_y2(const double) const;
            const double										m_y3(const double) const;
            const double                                        m_dy3(const double) const;
            void                                                compD(const Point&, const Point&, const Point&);
            void                                                compNormal(const Point&, const Point&, const Point&);
            const int											createS();
            //std::vector<double>									m_w[m_number];
            //std::vector<CFESolution>							m_w{ 4 };
            int													m_s;
            int													m_sp;
        };

		class CRectangleConstantBasis : public CShapeFunction<double>
		{
		public:
			CRectangleConstantBasis();
			CRectangleConstantBasis(const Point&, const Point&, const Point&, const Point&, const int order);
			CRectangleConstantBasis(const Point*, const int order);
			CRectangleConstantBasis(const CRectangleConstantBasis&);
			CRectangleConstantBasis& operator=(const CRectangleConstantBasis& t)
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
				return *this;
			}
			~CRectangleConstantBasis() {};
			const int				                            GetNumberOfShapeFunctions() const;
			//const DForm<0>*				                        GetShapeFunction(const int, const Point&) const;
			const double				                        GetShapeFunction(const int, const Point&) const;
			const Point											GetGradShapeFunction(const int, const Point&) const;
			const Point											GetNormal() const;
			void												ReverseNormal();
			const double										GetValue(const Point&) const;
			const int											IncreaseOrder();
			//const int											SetValue(const int, CSolution* value);
			//CSolution*										GetValue(const unsigned int);
			//const CFESolution									GetValue(const int) const;
			const double										GetMeasure() const { return m_det; };
			//const unsigned int									GetOrder() const;
			//const std::function<const DForm<0>*(const Point&)>  GetShapeFunction(const int) const;
		private:
			int													m_number;
			int													m_order;
			int													m_1dorder;
			Point                                               m_normal;
			std::vector<Mesh::Point>							m_points;
			double                                              m_det;
			double												m_hx, m_hy;
			const double										m_x1(const double) const;
			const double										m_x2(const double) const;
			const double										m_y1(const double) const;
			const double										m_y2(const double) const;
			void                                                compD(const Point&, const Point&, const Point&);
			void                                                compNormal(const Point&, const Point&, const Point&);
			const int											createS();
			//std::vector<double>									m_w[m_number];
			//std::vector<CFESolution>							m_w{ 1 };
			int													m_s;
			int													m_sp;
		};
	}
}
#endif // CORENC_MESH_RECTANGLE_H_
