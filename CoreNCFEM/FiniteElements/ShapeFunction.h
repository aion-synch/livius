#ifndef CORENC_MESH_ShapeFunction_h
#define CORENC_MESH_ShapeFunction_h
#include "../Point.h"
#include <functional>
#include "../FESolution.h"
namespace corenc
{
	namespace Mesh
	{
		template<class Type>
		class CShapeFunction
		{
		public:
			CShapeFunction() {}
			CShapeFunction(const Point*) {}
			virtual ~CShapeFunction() {}
			virtual const int				GetNumberOfShapeFunctions() const = 0;
			//virtual const  std::function<const DiffForm*(const Point&)> GetShapeFunction(const int) const = 0;
			//virtual const DiffForm*  GetShapeFunction(const int) const = 0;
			virtual const double			GetShapeFunction(const int, const Point&) const = 0;
			virtual const Point				GetGradShapeFunction(const int, const Point&) const = 0;
			virtual const Point				GetNormal() const = 0;
			virtual void					ReverseNormal() = 0;
			virtual const double			GetMeasure() const = 0;
			//virtual const Type				GetValue(const Point&) const = 0;
			//virtual const int				SetValue(const unsigned int, const Type& value) = 0;
			//virtual const Type				GetValue(const unsigned int) const = 0;
			//virtual const int				SetValue(const )
			//virtual CSolution*		GetValue(const unsigned int) = 0;
			//virtual const int				SetValue(const int, CSolution*) = 0;
		};
	}
}
#endif /* CORENC_MESH_ShapeFunction_h */
