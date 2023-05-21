#ifndef CORENC_MESH_Shape_h
#define CORENC_MESH_Shape_h
#include <functional>
#include <vector>
#include "../Point.h"
namespace corenc
{
	using scalar_func = std::function<const double(const Mesh::Point&)>;
	using vector_func = std::function<const Mesh::Point(const Mesh::Point&)>;
	namespace Mesh
	{
		//class Point;
		enum class NODES
		{
			FIRST,
			LAST
		};
		class CShape
		{
		public:
			CShape() {}
			CShape(const int*) {}
			virtual ~CShape() {}
			virtual const int					GetNumberOfNodes() const { return 0; };
			virtual const int					GetNumberOfEdges() const { return 0; };
			virtual const int					GetNumberOfFacets() const { return 0; };
			virtual const int     				GetNode(const int) const { return 1; };
			virtual const int     				GetNode(const NODES&) const { return 1; };
			virtual const int     				GetEdge(const int) const { return -1; };
			virtual const int     				GetFacet(const int) const { return -1; };
			virtual const double  				Integrate(const scalar_func&, const std::vector<Point>&) const = 0;
			virtual const Point  				Integrate(const vector_func&, const std::vector<Point>&) const = 0;
			virtual const std::vector<double>	Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const = 0;
			virtual void          				SetNode(const int, const int) = 0;
			virtual void          				SetEdge(const int, const int) {};
			virtual void          				SetFacet(const int, const int) {};
			//virtual std::istream&			operator>>(std::istream&) = 0;
		};
	}

}
#endif /* CORENC_MESH_Shape_h */
