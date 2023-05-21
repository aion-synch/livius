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
		class CTriangularMeshLinear : public CMesh<>
		{
		public:
			CTriangularMeshLinear();
			CTriangularMeshLinear(const std::string& file_name);
			CTriangularMeshLinear(const CTriangularMeshLinear&);
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
			~CTriangularMeshLinear();
		private:
			std::vector<CElement<>*>			m_elems;
			std::vector<CElement<>*>			m_edges;
			std::vector<Point>					m_points;
			std::map<int, CParameter>			m_params;
		public:
			auto								GetElements() -> decltype(m_elems) { return m_elems; };
			auto								GetBoundary() -> decltype(m_edges) { return m_edges; };
		};
	}
}
#endif /* CORENC_MESH_TriangularMesh_h */

