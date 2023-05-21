#ifndef CORENC_MESH_Mesh_h
#define CORENC_MESH_Mesh_h
#include "FiniteElements/FiniteElement.h"
namespace corenc
{
	namespace Mesh
	{
		enum Meshes
		{
			Mesh1D = 0,
			TriangularMesh = 1,
			TetrahedralMesh = 2
		};
		template<class T = bool>
		class CMesh;
		template<class T>
		class CMesh
		{
		public:
			CMesh() {}
			virtual ~CMesh() {}
			virtual const unsigned int				GetNumberOfNodes() const = 0;
			virtual const unsigned int				GetNumberOfElements() const = 0;
			virtual const int						FindElement(const Point&) const = 0;
			virtual const unsigned int				GetNumberOfBoundaries() const = 0;
			virtual const CElement<T>*				GetElement(const unsigned int) const = 0;
			virtual const CElement<T>*				GetBoundary(const unsigned int) const = 0;
			virtual const Point						GetNode(const unsigned int) const = 0;
			virtual const double					getSolution(const unsigned int element, const unsigned int node) const = 0;
			virtual const int						updateSolution(const unsigned int element, const unsigned int node, const double value) = 0;
			virtual const std::vector<double>		getSolution() const = 0;
			virtual const int						updateSolution(const std::vector<double>&) = 0;
			//virtual const int						updateSolution(const std::vector<CSolution*>&) = 0;
			virtual const int						updateSolution(const unsigned int element, const unsigned int node, CSolution* value) = 0;
			virtual const double					getParameter(Parameters, const unsigned int, const Point&) const = 0;
			virtual const double					getParameter(Parameters, const unsigned int, const int) const = 0;
			virtual const int						setParameter(Parameters, const double, const unsigned int) = 0;
			virtual const double					getMinSize() const = 0;
			virtual const int						updateSolution(const unsigned int node, const double value) = 0;
		};
		template<>
		class CMesh<bool>
		{
		public:
			CMesh() {}
			virtual ~CMesh() {}
			virtual const unsigned int				GetNumberOfNodes() const = 0;
			virtual const unsigned int				GetNumberOfElements() const = 0;
			virtual const int						FindElement(const Point&) const = 0;
			virtual const unsigned int				GetNumberOfBoundaries() const = 0;
			virtual const CElement<>*				GetElement(const unsigned int) const = 0;
			virtual const CElement<>*				GetBoundary(const unsigned int) const = 0;
			virtual const Point						GetNode(const unsigned int) const = 0;
			virtual const double					getSolution(const unsigned int element, const unsigned int node) const = 0;
			virtual const int						updateSolution(const unsigned int element, const unsigned int node, const double value) = 0;
			virtual const std::vector<double>		getSolution() const = 0;
			virtual const int						updateSolution(const std::vector<double>&) = 0;
			//virtual const int						updateSolution(const std::vector<CSolution*>&) = 0;
			virtual const int						updateSolution(const unsigned int element, const unsigned int node, CSolution* value) = 0;
			virtual const double					getParameter(Parameters, const unsigned int, const Point& p) const = 0;
			virtual const double					getParameter(Parameters, const unsigned int, const int) const = 0;
			virtual const int						setParameter(Parameters, const double, const unsigned int) = 0;
			virtual const double					getMinSize() const = 0;
		};
	}
}

#endif /* CORENC_MESH_Mesh_h */
