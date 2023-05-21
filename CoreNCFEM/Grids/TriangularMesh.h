#ifndef CORENC_MESH_TriangularMesh_h
#define CORENC_MESH_TriangularMesh_h
#include "../Mesh.h"

#include "../Point.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "../Parameter.h"
namespace corenc
{
	namespace Mesh
	{
		class CTriangularMesh : public CMesh<>
		{
		public:
			CTriangularMesh();
			CTriangularMesh(const std::string& file_name);
			CTriangularMesh(const CTriangularMesh&);
			CTriangularMesh(const Point& p1, const Point& p2, const int nx, const int ny);
			CTriangularMesh& operator=(const CTriangularMesh& tr)
			{
				const int sz_el = tr.m_elems.size();
				const int sz_pt = tr.m_points.size();
				const int sz_bpt = tr.m_basepoints.size();
				const int sz_bel = tr.m_elemsbase.size();
				const int sz_ed = tr.m_edges.size();
				const int sz_bed = tr.m_edgesbase.size();
				m_elems.resize(sz_el);
				m_edges.resize(sz_ed);
				m_points.resize(sz_pt);
				m_basepoints.resize(sz_bpt);
				m_elemsbase.resize(sz_bel);
				m_edgesbase.resize(sz_bed);
				int i = 0;
				for (i = 0; i < sz_el; ++i)
					m_elems[i] = tr.m_elems[i]->Clone();
				for (i = 0; i < sz_ed; ++i)
					m_edges[i] = tr.m_edges[i]->Clone();
				for (i = 0; i < sz_pt; ++i)
					m_points[i] = tr.m_points[i];
				for (i = 0; i < sz_bpt; ++i)
					m_basepoints[i] = tr.m_basepoints[i];
				for (i = 0; i < sz_bel; ++i)
					m_elemsbase[i] = tr.m_elemsbase[i]->Clone();
				for (i = 0; i < sz_bed; ++i)
					m_edgesbase[i] = tr.m_edgesbase[i]->Clone();
				m_params = tr.m_params;
                m_bnds.resize(tr.m_bnds.size());
                for (i = 0; i < m_bnds.size(); ++i)
                    m_bnds[i] = tr.m_bnds[i];
				m_order = tr.m_order;
                m_offsets = tr.m_offsets;
				return *this;
			}
			CTriangularMesh*			Clone() const
			{
				return new CTriangularMesh(*this);
			};
			//const bool operator<(const CTriangularMesh& mesh) const
			//{
			//	if(m_points.size() < mesh.m_points.size())
			//		return true;
			//	return false;
			//}
			const unsigned int					GetNumberOfElements() const;
			const unsigned int					GetNumberOfNodes() const;
			const unsigned int					GetNumberOfBoundaries() const;
			const int							FindElement(const Point&) const;
			const Point							GetNode(const unsigned int) const;
			const CElement<>* 					GetElement(const unsigned int) const;
			const CElement<>*					GetBoundary(const unsigned int) const;
			const double						getMinSize() const { return 0.; };
			const double						getSolution(const unsigned int element, const unsigned int node) const;
			const int							updateSolution(const unsigned int element, const unsigned int node, const double value);
			const std::vector<double>			getSolution() const;
			const int							updateSolution(const std::vector<double>&);
			const int							updateSolution(const unsigned int element, const unsigned int node, CSolution* value);
			const double						getParameter(Parameters, const unsigned int, const Point&) const;
			const double						getParameter(Parameters, const unsigned int, const int) const;
			const int							setParameter(Parameters, const double, const unsigned int);
			const int							setParameter(const CParameter&, const unsigned int type);
			const int							updateSolution(const unsigned int node, const double value);
			const int							refine_h();
			const int							refine_p();
			const int							refine_hp();
            const int                           set4thOrder();
            const int                           set2ndOrder();
            const int                           set3rdOrder();
            const int                           interpolate(const int node) const;
            const int                           GetNumberOfINodes() const;
			~CTriangularMesh();
		private:
			std::vector<CElement<>*>			m_elems;
			std::vector<CElement<>*>			m_edges;
			std::vector<Point>					m_points;
			std::vector<Point>					m_basepoints;
			std::vector<CElement<>*>			m_elemsbase;
			std::vector<CElement<>*>			m_edgesbase;
			std::map<int, CParameter>			m_params;
            std::vector<int>                    m_offsets;
            std::vector<int>                    m_bnds;
			const double						CompSquare(const Point& p1, const Point& p2, const Point& p3) const;
            void                                set3rdNodes();
            void                                set4thNodes_1();
            void                                set4thNodes_2();
			int									m_order;
		public:
			auto								GetElements() -> decltype(m_elems) { return m_elems; };
			auto								GetBoundary() -> decltype(m_edges) { return m_edges; };
		};
	}
}
#endif /* CORENC_MESH_TriangularMesh_h */

