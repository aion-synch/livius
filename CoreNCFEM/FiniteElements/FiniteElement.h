#ifndef CORENC_MESH_FINITEELEMENT_H_
#define CORENC_MESH_FINITEELEMENT_H_

#include <functional>
#include <iostream>
#include <vector>
#include "../Point.h"
#include "../FESolution.h"
namespace corenc
{
	namespace Mesh
	{
		using function_dp = std::function<const double(const Point&)>;
		enum Elements
		{
			Interval = 0,
			Triangle = 1,
			Rectangle = 2,
			Tetrahedron = 3,
			Cube = 4
		};
		
		template<class T = bool>
		class CElement;

		template<>
		class CElement<bool>
		{
		public:
			CElement() {}
			virtual ~CElement() {}
			virtual const int									GetType() const = 0;
			virtual CElement<>*									Clone() const = 0;
			virtual const int									GetDoFs() const = 0;
			virtual const int									GetNode(const int) const = 0;
			virtual const int									GetNeighbour(const int) const = 0;
			virtual void										SetNeighbour(const int k, const int elem) = 0;
			virtual void										SetType(const int) = 0;
			virtual void										SetNode(const int, const int) = 0;
			virtual const int									GetNumberOfNodes() const = 0;
			//virtual void            SetShapeFunction(const unsigned int, const std::function<const DiffForm(const Point&)>&) = 0;
			//virtual const DiffForm* GetShapeFunction(const int, const Point&) const = 0;
			//virtual const double								GetShapeFunction(const unsigned int, const std::vector<double>&) const = 0;
			virtual const double								GetShapeFunction(const int, const Point&) const = 0;
			virtual const Point									GetGradShapeFunction(const int, const Point&) const = 0;
			//virtual const std::vector<double>					GetGradShapeFunction(const unsigned int, const std::vector<double>&) const = 0;
			virtual const Point									GetNormal() const = 0;
			//virtual const std::vector<double>					GetNormal() const = 0;
			virtual void										ReverseNormal() = 0;
			virtual const double								Integrate(const function_dp&, const std::vector<Point>& v) const = 0;
			//virtual const double								Integrate(const std::function<const double(const std::vector<double>&)>&, const std::vector<std::vector<double>>& v) const = 0;
			virtual const Point									Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const = 0;
			//virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const std::vector<double>&)>&, const std::vector<std::vector<double>>& v) const = 0;
			virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const = 0;
			virtual const double								GetWeight(const int, const std::vector<Point>& verts, const function_dp& f) const = 0;
			//virtual const Type									GetValue(const unsigned number) const = 0;
			//virtual const Type									GetValue(const Mesh::Point&) const = 0;
			//virtual const int									SetValue(const unsigned int number, const Type& value) = 0;
			//virtual const int									SetValue(const int number, CSolution* value) = 0;
			//virtual CSolution*									GetValue(const int) = 0;
			virtual const int									IncreaseOrder() = 0;
			virtual const double								GetMeasure() const = 0;
			//virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const std::vector<double>&)>&, const std::vector<std::vector<double>>&) const = 0;
			//virtual std::function<const DiffForm(const Point&)> GetShapeFunction(const int) = 0;
			//virtual const double	GetShapeFunction(const int, const Point&) const = 0;
		};
		template<class T>
		class CElement
		{
		public:
			CElement() {}
			virtual ~CElement() {}
			virtual const int									GetType() const = 0;
			virtual CElement*									Clone() const = 0;
			virtual const int									GetDoFs() const = 0;
			virtual const int									GetNode(const int) const = 0;
			virtual const int									GetNeighbour(const int) const = 0;
			virtual void										SetNeighbour(const int k, const int elem) = 0;
			virtual void										SetType(const int) = 0;
			virtual void										SetNode(const int, const int) = 0;
			virtual const int									GetNumberOfNodes() const = 0;
			//virtual void            SetShapeFunction(const unsigned int, const std::function<const DiffForm(const Point&)>&) = 0;
			//virtual const DiffForm* GetShapeFunction(const int, const Point&) const = 0;
			//virtual const double								GetShapeFunction(const unsigned int, const std::vector<double>&) const = 0;
			virtual const double								GetShapeFunction(const int, const Point&) const = 0;
			virtual const Point									GetGradShapeFunction(const int, const Point&) const = 0;
			//virtual const std::vector<double>					GetGradShapeFunction(const unsigned int, const std::vector<double>&) const = 0;
			virtual const Point									GetNormal() const = 0;
			//virtual const std::vector<double>					GetNormal() const = 0;
			virtual void										ReverseNormal() = 0;
			virtual const int									IncreaseOrder() = 0;
			virtual const double								Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const = 0;
			//virtual const double								Integrate(const std::function<const double(const std::vector<double>&)>&, const std::vector<std::vector<double>>& v) const = 0;
			virtual const Point									Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const = 0;
			//virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const std::vector<double>&)>&, const std::vector<std::vector<double>>& v) const = 0;
			virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const = 0;
			//virtual const Type									GetValue(const unsigned number) const = 0;
			//virtual const Type									GetValue(const Mesh::Point&) const = 0;
			//virtual const int									SetValue(const unsigned int number, const Type& value) = 0;
			//virtual const int									SetValue(const unsigned int number, CSolution* value) = 0;
			//virtual CSolution*							GetValue(const unsigned int) = 0;
			virtual const double								GetMeasure() const = 0;
			virtual const double								GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const = 0;
			//virtual const T										GetValue(const int number) const = 0;
			//virtual const T										GetValue(const Point& p) const = 0;
			//virtual const int									SetValue(const int number, const T& value) = 0;
			//virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const std::vector<double>&)>&, const std::vector<std::vector<double>>&) const = 0;
			//virtual std::function<const DiffForm(const Point&)> GetShapeFunction(const int) = 0;
			//virtual const double	GetShapeFunction(const int, const Point&) const = 0;
		};




		// Set of nodes
		// Set of shape function
		// Set of degrees of freedom ; don't use pls
		// Type of the weights aligned with the degrees of freedom
		// The weights should be inside of the set of shape functions and the types should be same
		template<class Shape, class ShapeFunction, class DoF = bool, class T = bool>
		class CFiniteElement;
		template<class Shape, class ShapeFunction, class DoF, class T>
		class CFiniteElement: public CElement<T>
		{
		public:
			CFiniteElement() {}
			CFiniteElement(const int* nodes, const Point* points, const int dofs) :
				m_shape{ nodes },
				m_shapefunctions{ points },
				m_dofs{ dofs },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const int* nodes, const Point* points) :
				m_shape{ nodes },
				m_shapefunctions{ points },
				m_dofs{ 0 },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& f, const DoF& d) :
				m_shape{ shape },
				m_shapefunctions{ f },
				m_dofs{ d },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& shfunc, const DoF& dofs, const int type) :
				m_shape{ shape },
				m_shapefunctions{ shfunc },
				m_dofs{ dofs },
				m_type{ type } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
                        CFiniteElement(const Shape& shape, const ShapeFunction& shfunc, const DoF& dofs, const int type, const int* neigs) :
                                m_shape{ shape },
                                m_shapefunctions{ shfunc },
                                m_dofs{ dofs },
                                m_type{ type } {
                                m_neighbours[0] = neigs[0]; m_neighbours[1] = neigs[1];
                        };
			CFiniteElement(const CFiniteElement<Shape, ShapeFunction, DoF>& e) :
				m_shape{ e.m_shape },
				m_shapefunctions{ e.m_shapefunctions },
				m_dofs{ e.m_dofs },
				m_type{ e.m_type } {
				m_neighbours[0] = e.m_neighbours[0]; m_neighbours[1] = e.m_neighbours[1];
			};
			CElement<T>*			Clone() const
			{
                                return new CFiniteElement<Shape, ShapeFunction, DoF, T>(m_shape, m_shapefunctions, m_dofs, m_type, m_neighbours);
			};
			friend const bool		operator==(const CFiniteElement& e1, const CFiniteElement& e2)
			{
				if (e1.m_shape == e2.m_shape)
					return true;
				return false;
			}
			~CFiniteElement() {}
			const int               GetType() const;
			const int               GetNode(const int) const;
			const int				GetNeighbour(const int) const;
			const Shape             GetShape() const;
			const ShapeFunction     GetShapeFunctions() const;
			const DoF               GetDoF() const;
			const int				GetDoFs() const;
			void					SetNeighbour(const int k, const int elem);
			void                    SetType(const int);
			void                    SetShapeFunction(const int, const ShapeFunction&);
			void                    SetDoF(const DoF&);
			void                    SetShape(const Shape&);
			const int				IncreaseOrder();
			void                    SetNode(const int, const int);
			const int				GetNumberOfNodes() const;
			//const int				SetValue(const int number, CSolution* value);
			const double			GetMeasure() const;
			//CSolution*				GetValue(const int);
			//void                    SetShapeFunction(const int, const std::function<const DiffForm(const Point&)>&);
			const double			Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point				Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const std::vector<double>	Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
			//const std::function<const DiffForm(const Point&)>          GetShapeFunction(const int) const;
			//const DiffForm*			GetShapeFunction(const int, const Point&);
			const double			GetShapeFunction(const int, const Point&) const;
			const Point				GetGradShapeFunction(const int, const Point&) const;
			const Point				GetNormal() const;
			void					ReverseNormal();
			const double			GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
			//const T					GetValue(const int number) const;
			//const T					GetValue(const Point& p) const;
			//const int				SetValue(const int number, const T& value);
			CFiniteElement&         operator=(const CFiniteElement& e)
			{
				m_shape = e.m_shape;
				m_shapefunctions = e.m_shapefunctions;
				m_dofs = e.m_dofs;
				m_type = e.m_type;
				return *this;
			}
			friend std::istream&	operator>>(std::istream& is, CFiniteElement& k)
			{
				is >> k.m_shape;
				return is;
			}
			//const DiffForm          GetDShapeFunction(const int, const Point&);
		private:
			Shape                   m_shape;
			ShapeFunction           m_shapefunctions;
			DoF                     m_dofs;
			int                     m_type;
			int						m_neighbours[2];
		};

		template<class Shape, class ShapeFunction, class DoF>
		class CFiniteElement<Shape, ShapeFunction, DoF, bool> : public CElement<>
		{
		public:
			CFiniteElement() {}
			CFiniteElement(const int* nodes, const Point* points, const int dofs) :
				m_shape{ nodes },
				m_shapefunctions{ points },
				m_dofs{ dofs },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const int* nodes, const Point* points) :
				m_shape{ nodes },
				m_shapefunctions{ points },
				m_dofs{ 0 },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& f, const DoF& d) :
				m_shape{ shape },
				m_shapefunctions{ f },
				m_dofs{ d },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& shfunc, const DoF& dofs, const int type) :
				m_shape{ shape },
				m_shapefunctions{ shfunc },
				m_dofs{ dofs },
				m_type{ type } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& shfunc, const DoF& dofs, const int type, const int* neigh) :
				m_shape{ shape },
				m_shapefunctions{ shfunc },
				m_dofs{ dofs },
				m_type{ type } {
				m_neighbours[0] = neigh[0]; m_neighbours[1] = neigh[1];
			};
			CFiniteElement(const CFiniteElement<Shape, ShapeFunction, DoF>& e) :
				m_shape{ e.m_shape },
				m_shapefunctions{ e.m_shapefunctions },
				m_dofs{ e.m_dofs },
				m_type{ e.m_type } {
				m_neighbours[0] = e.m_neighbours[0]; m_neighbours[1] = e.m_neighbours[1];
			};
			friend const bool		operator==(const CFiniteElement& e1, const CFiniteElement& e2)
			{
				if (e1.m_shape == e2.m_shape)
					return true;
				return false;
			}
			// don't forget to delete after the call
			CElement<>*			Clone() const
			{
				return new CFiniteElement<Shape, ShapeFunction, DoF>(m_shape, m_shapefunctions, m_dofs, m_type, m_neighbours);
			};
			~CFiniteElement() {}
			const int               GetType() const;
			const int               GetNode(const int) const;
			const int				GetNeighbour(const int) const;
			const Shape             GetShape() const;
			const ShapeFunction     GetShapeFunctions() const;
			const DoF               GetDoF() const;
			const int				GetDoFs() const;
			void					SetNeighbour(const int k, const int elem);
			void                    SetType(const int);
			void                    SetShapeFunction(const int, const ShapeFunction&);
			void                    SetDoF(const DoF&);
			void                    SetShape(const Shape&);
			void                    SetNode(const int, const int);
			const int				GetNumberOfNodes() const;
			//const int				SetValue(const int number, CSolution* value);
			const int				IncreaseOrder();
			const double			GetMeasure() const;
			//CSolution*				GetValue(const int);
			//void                    SetShapeFunction(const int, const std::function<const DiffForm(const Point&)>&);
			const double			Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point				Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const std::vector<double>	Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
			//const std::function<const DiffForm(const Point&)>          GetShapeFunction(const int) const;
			//const DiffForm*			GetShapeFunction(const int, const Point&);
			const double			GetShapeFunction(const int, const Point&) const;
			const Point				GetGradShapeFunction(const int, const Point&) const;
			const Point				GetNormal() const;
			void					ReverseNormal();
			const double			GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
			CFiniteElement&         operator=(const CFiniteElement& e)
			{
				m_shape = e.m_shape;
				m_shapefunctions = e.m_shapefunctions;
				m_dofs = e.m_dofs;
				m_type = e.m_type;
				return *this;
			}
			friend std::istream&	operator>>(std::istream& is, CFiniteElement& k)
			{
				is >> k.m_shape;
				return is;
			}
			//const DiffForm          GetDShapeFunction(const int, const Point&);
		private:
			Shape                   m_shape;
			ShapeFunction           m_shapefunctions;
			DoF                     m_dofs;
			int                     m_type;
			int						m_neighbours[2];
		};


		template<class Shape, class ShapeFunction>
		class CFiniteElement<Shape, ShapeFunction, bool, bool> : public CElement<>
		{
		public:
			CFiniteElement() {}
			CFiniteElement(const int* nodes, const Point* points, const int dofs) :
				m_shape{ nodes },
				m_shapefunctions{ points, dofs },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const int* nodes, const Point* points, const int dofs, const int type) :
				m_shape{ nodes, dofs },
				m_shapefunctions{ points, dofs },
				m_type{ type } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const int* nodes, const Point* points) :
				m_shape{ nodes },
				m_shapefunctions{ points },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& f) :
				m_shape{ shape },
				m_shapefunctions{ f },
				m_type{ -1 } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
			CFiniteElement(const Shape& shape, const ShapeFunction& shfunc, const int type) :
				m_shape{ shape },
				m_shapefunctions{ shfunc },
				m_type{ type } {
				m_neighbours[0] = -1; m_neighbours[1] = -1;
			};
                        CFiniteElement(const Shape& shape, const ShapeFunction& shfunc, const int type, const int* neigs) :
                                m_shape{ shape },
                                m_shapefunctions{ shfunc },
                                m_type{ type } {
                                m_neighbours[0] = neigs[0]; m_neighbours[1] = neigs[1];
                        };
			CFiniteElement(const CFiniteElement<Shape, ShapeFunction>&e) :
				m_shape{ e.m_shape },
				m_shapefunctions{ e.m_shapefunctions },
				m_type{ e.m_type } {
				m_neighbours[0] = e.m_neighbours[0]; m_neighbours[1] = e.m_neighbours[1];
			};
			friend const bool		operator==(const CFiniteElement& e1, const CFiniteElement& e2)
			{
				if (e1.m_shape == e2.m_shape)
					return true;
				return false;
			}
			// don't forget to delete after the call
			CElement<>*			Clone() const
			{
                                return new CFiniteElement<Shape, ShapeFunction>(m_shape, m_shapefunctions, m_type, m_neighbours);
			};
			~CFiniteElement() {}
			const int               GetType() const;
			const int               GetNode(const int) const;
			const int				GetNeighbour(const int) const;
			const Shape             GetShape() const;
			const ShapeFunction     GetShapeFunctions() const;
			const int				GetDoFs() const;
			void					SetNeighbour(const int k, const int elem);
			void                    SetType(const int);
			void                    SetShapeFunction(const int, const ShapeFunction&);
			void                    SetShape(const Shape&);
			void                    SetNode(const int, const int);
			const int				GetNumberOfNodes() const;
			//const int				SetValue(const int number, CSolution* value);
			const int				IncreaseOrder();
			const double			GetMeasure() const;
			//CSolution*				GetValue(const int);
			//void                    SetShapeFunction(const int, const std::function<const DiffForm(const Point&)>&);
			const double			Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const;
			const Point				Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const;
			const std::vector<double>	Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const;
			//const std::function<const DiffForm(const Point&)>          GetShapeFunction(const int) const;
			//const DiffForm*			GetShapeFunction(const int, const Point&);
			const double			GetShapeFunction(const int, const Point&) const;
			const Point				GetGradShapeFunction(const int, const Point&) const;
			const Point				GetNormal() const;
			void					ReverseNormal();
			const double			GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const;
			CFiniteElement&         operator=(const CFiniteElement& e)
			{
				m_shape = e.m_shape;
				m_shapefunctions = e.m_shapefunctions;
				m_type = e.m_type;
				return *this;
			}
			friend std::istream&	operator>>(std::istream& is, CFiniteElement& k)
			{
				is >> k.m_shape;
				return is;
			}
			//const DiffForm          GetDShapeFunction(const int, const Point&);
		private:
			Shape                   m_shape;
			ShapeFunction           m_shapefunctions;
			int                     m_type;
			int						m_neighbours[2];
		};





		// implementation template<class Shape, class ShapeFunction, class DoF>
		// CFiniteElement<Shape, ShapeFunction, DoF>
		template<class Shape, class ShapeFunction, class DoF>
		const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetType() const
		{
			return m_type;
		}
		//template<class Shape, class ShapeFunction, class DoF>
		//const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetValue(const int number, CSolution* value)
		//{
		//	return m_shapefunctions.SetValue(number, value);
		//}
		//template<class Shape, class ShapeFunction, class DoF>
		//CSolution* CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetValue(const int p)
		//{
		//	auto&& val = m_shapefunctions.GetValue(p);
		//	return const_cast<CSolution*>(static_cast<const CSolution*>(&val));;
		//}
		template<class Shape, class ShapeFunction, class DoF>
		const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetNode(const int k) const
		{
			return m_shape.GetNode(k);
		}
		template<class Shape, class ShapeFunction, class DoF>
		const double CFiniteElement<Shape, ShapeFunction, DoF, bool>::Integrate(const std::function<const double(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}

		template<class Shape, class ShapeFunction, class DoF>
		const Point CFiniteElement<Shape, ShapeFunction, DoF, bool>::Integrate(const std::function<const Point(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}

		template<class Shape, class ShapeFunction, class DoF>
		const std::vector<double> CFiniteElement<Shape, ShapeFunction, DoF, bool>::Integrate(const std::function<const std::vector<double>(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}
		template<class Shape, class ShapeFunction, class DoF>
		const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetNeighbour(const int k) const
		{
			return m_neighbours[k];
		}

		template<class Shape, class ShapeFunction, class DoF>
		const Point CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetNormal() const
		{
			return m_shapefunctions.GetNormal();
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::ReverseNormal()
		{
			m_shapefunctions.ReverseNormal();
		}

		template<class Shape, class ShapeFunction, class DoF>
		const double CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetMeasure() const
		{
			return m_shapefunctions.GetMeasure();
		}

		//template<class Shape, class ShapeFunction, class DoF, class T>
		//inline const T CFiniteElement<Shape, ShapeFunction, DoF, T>::GetValue(const int number) const
		//{
		//	return m_shapefunctions.GetValue(number);
		//}

		//template<class Shape, class ShapeFunction, class DoF, class T>
		//inline const T CFiniteElement<Shape, ShapeFunction, DoF, T>::GetValue(const Point & p) const
		//{
		//	return m_shapefunctions.GetValue(p);
		//}

		//template<class Shape, class ShapeFunction, class DoF, class T>
		//inline const int CFiniteElement<Shape, ShapeFunction, DoF, T>::SetValue(const int number, const T & value)
		//{
		//	return m_shapefunctions.SetValue(number, value);
		//}

		template<class Shape, class ShapeFunction, class DoF>
		const Shape CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetShape() const
		{
			return m_shape;
		}

		template<class Shape, class ShapeFunction, class DoF>
		const ShapeFunction CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetShapeFunctions() const
		{
			return m_shapefunctions;
		}

		template<class Shape, class ShapeFunction, class DoF>
		const DoF CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetDoF() const
		{
			return m_dofs;
		}

		template<class Shape, class ShapeFunction, class DoF>
		const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetNumberOfNodes() const
		{
			return m_shape.GetNumberOfNodes();
		}

		template<class Shape, class ShapeFunction, class DoF>
		const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetDoFs() const
		{
			return m_shapefunctions.GetNumberOfShapeFunctions();
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetNeighbour(const int k, const int elem)
		{
			m_neighbours[k] = elem;
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetShapeFunction(const int k, const ShapeFunction& func)
		{
			m_shapefunctions = func;
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetDoF(const DoF& dof)
		{
			m_dofs = dof;
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetShape(const Shape &shape)
		{
			m_shape = shape;
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetType(const int k)
		{
			m_type = k;
		}

		template<class Shape, class ShapeFunction, class DoF>
		void CFiniteElement<Shape, ShapeFunction, DoF, bool>::SetNode(const int k, const int node)
		{
			m_shape.SetNode(k, node);
		}
		template<class Shape, class ShapeFunction, class DoF>
		const double CFiniteElement<Shape, ShapeFunction, DoF>::GetShapeFunction(const int k, const Mesh::Point &p) const
		{
			return m_shapefunctions.GetShapeFunction(k, p);
		}

		template<class Shape, class ShapeFunction, class DoF>
		const Point CFiniteElement<Shape, ShapeFunction, DoF, bool>::GetGradShapeFunction(const int k, const Mesh::Point &p) const
		{
			return m_shapefunctions.GetGradShapeFunction(k, p);
		}

		template<class Shape, class ShapeFunction, class DoF>
		const int CFiniteElement<Shape, ShapeFunction, DoF, bool>::IncreaseOrder()
		{
			if (m_shape.IncreaseOrder())
				return 1;
			if (m_shapefunctions.IncreaseOrder())
				return 1;
			return 0;
		}
		template<class Shape, class ShapeFunction, class DoF>
		const double CFiniteElement<Shape, ShapeFunction, DoF>::GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
		{
			return 0.0;
		}
		// fin.

		// implementation template<class Shape, class ShapeFunction, class DoF, class T>
		// CFiniteElement<Shape, ShapeFunction, DoF, T>
		template<class Shape, class ShapeFunction, class DoF, class T>
		const int CFiniteElement<Shape, ShapeFunction, DoF, T>::GetType() const
		{
			return m_type;
		}
		template<class Shape, class ShapeFunction, class DoF, class T>
		const double CFiniteElement<Shape, ShapeFunction, DoF, T>::GetMeasure() const
		{
			return m_shapefunctions.GetMeasure();
		}
		//template<class Shape, class ShapeFunction, class DoF, class T>
		//const int CFiniteElement<Shape, ShapeFunction, DoF, T>::SetValue(const int number, CSolution* value)
		//{
		//	return m_shapefunctions.SetValue(number, value);
		//}
		//template<class Shape, class ShapeFunction, class DoF, class T>
		//CSolution* CFiniteElement<Shape, ShapeFunction, DoF, T>::GetValue(const int p)
		//{
		//	auto&& val = m_shapefunctions.GetValue(p);
		//	return const_cast<CSolution*>(static_cast<const CSolution*>(&val));
		//}
		template<class Shape, class ShapeFunction, class DoF, class T>
		const int CFiniteElement<Shape, ShapeFunction, DoF, T>::GetNode(const int k) const
		{
			return m_shape.GetNode(k);
		}
		template<class Shape, class ShapeFunction, class DoF, class T>
		const double CFiniteElement<Shape, ShapeFunction, DoF, T>::Integrate(const std::function<const double(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const Point CFiniteElement<Shape, ShapeFunction, DoF, T>::Integrate(const std::function<const Point(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const std::vector<double> CFiniteElement<Shape, ShapeFunction, DoF, T>::Integrate(const std::function<const std::vector<double>(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}
		template<class Shape, class ShapeFunction, class DoF, class T>
		const int CFiniteElement<Shape, ShapeFunction, DoF, T>::GetNeighbour(const int k) const
		{
			return m_neighbours[k];
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const Point CFiniteElement<Shape, ShapeFunction, DoF, T>::GetNormal() const
		{
			return m_shapefunctions.GetNormal();
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::ReverseNormal()
		{
			m_shapefunctions.ReverseNormal();
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
                const double CFiniteElement<Shape, ShapeFunction, DoF, T>::GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
		{
                        return m_shapefunctions.GetWeight(node, verts, f);
		}


		template<class Shape, class ShapeFunction, class DoF, class T>
		const Shape CFiniteElement<Shape, ShapeFunction, DoF, T>::GetShape() const
		{
			return m_shape;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const ShapeFunction CFiniteElement<Shape, ShapeFunction, DoF, T>::GetShapeFunctions() const
		{
			return m_shapefunctions;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const DoF CFiniteElement<Shape, ShapeFunction, DoF, T>::GetDoF() const
		{
			return m_dofs;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const int CFiniteElement<Shape, ShapeFunction, DoF, T>::GetNumberOfNodes() const
		{
			return m_shape.GetNumberOfNodes();
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const int CFiniteElement<Shape, ShapeFunction, DoF, T>::GetDoFs() const
		{
			return m_shapefunctions.GetNumberOfShapeFunctions();
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::SetNeighbour(const int k, const int elem)
		{
			m_neighbours[k] = elem;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::SetShapeFunction(const int k, const ShapeFunction& func)
		{
			m_shapefunctions = func;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::SetDoF(const DoF& dof)
		{
			m_dofs = dof;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::SetShape(const Shape &shape)
		{
			m_shape = shape;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::SetType(const int k)
		{
			m_type = k;
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		void CFiniteElement<Shape, ShapeFunction, DoF, T>::SetNode(const int k, const int node)
		{
			m_shape.SetNode(k, node);
		}
		template<class Shape, class ShapeFunction, class DoF, class T>
		const double CFiniteElement<Shape, ShapeFunction, DoF, T>::GetShapeFunction(const int k, const Mesh::Point &p) const
		{
			return m_shapefunctions.GetShapeFunction(k, p);
		}

		template<class Shape, class ShapeFunction, class DoF, class T>
		const Point CFiniteElement<Shape, ShapeFunction, DoF, T>::GetGradShapeFunction(const int k, const Mesh::Point &p) const
		{
			return m_shapefunctions.GetGradShapeFunction(k, p);
		}
		template<class Shape, class ShapeFunction, class DoF, class T>
		const int CFiniteElement<Shape, ShapeFunction, DoF, T>::IncreaseOrder()
		{
			if (m_shape.IncreaseOrder())
				return 1;
			if (m_shapefunctions.IncreaseOrder())
				return 1;
			return 0;
		}

		// implementation template<class Shape, class ShapeFunction>
		// CFiniteElement<Shape, ShapeFunction>

		template<class Shape, class ShapeFunction>
		const int CFiniteElement<Shape, ShapeFunction, bool, bool>::GetType() const
		{
			return m_type;
		}
		template<class Shape, class ShapeFunction>
		const double CFiniteElement<Shape, ShapeFunction, bool, bool>::GetMeasure() const
		{
			return m_shapefunctions.GetMeasure();
		}
		//template<class Shape, class ShapeFunction>
		//const int CFiniteElement<Shape, ShapeFunction, bool, bool>::SetValue(const int number, CSolution* value)
		//{
		//	return m_shapefunctions.SetValue(number, value);
		//}
		//template<class Shape, class ShapeFunction>
		//CSolution* CFiniteElement<Shape, ShapeFunction, bool, bool>::GetValue(const int p)
		//{
		//	auto&& val = m_shapefunctions.GetValue(p);
		//	return const_cast<CSolution*>(static_cast<const CSolution*>(&val));;
		//}
		template<class Shape, class ShapeFunction>
		const int CFiniteElement<Shape, ShapeFunction, bool, bool>::GetNode(const int k) const
		{
			return m_shape.GetNode(k);
		}
		template<class Shape, class ShapeFunction>
		const double CFiniteElement<Shape, ShapeFunction, bool, bool>::Integrate(const std::function<const double(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}

		template<class Shape, class ShapeFunction>
		const Point CFiniteElement<Shape, ShapeFunction, bool, bool>::Integrate(const std::function<const Point(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}

		template<class Shape, class ShapeFunction>
		const std::vector<double> CFiniteElement<Shape, ShapeFunction, bool, bool>::Integrate(const std::function<const std::vector<double>(const Point&)>&f, const std::vector<Point>& v) const
		{
			return m_shape.Integrate(f, v);
		}
		template<class Shape, class ShapeFunction>
		const int CFiniteElement<Shape, ShapeFunction, bool, bool>::GetNeighbour(const int k) const
		{
			return m_neighbours[k];
		}

		template<class Shape, class ShapeFunction>
		const Point CFiniteElement<Shape, ShapeFunction, bool, bool>::GetNormal() const
		{
			return m_shapefunctions.GetNormal();
		}

		template<class Shape, class ShapeFunction>
		void CFiniteElement<Shape, ShapeFunction, bool, bool>::ReverseNormal()
		{
			m_shapefunctions.ReverseNormal();
		}


		template<class Shape, class ShapeFunction>
		const Shape CFiniteElement<Shape, ShapeFunction, bool, bool>::GetShape() const
		{
			return m_shape;
		}

		template<class Shape, class ShapeFunction>
		const ShapeFunction CFiniteElement<Shape, ShapeFunction, bool, bool>::GetShapeFunctions() const
		{
			return m_shapefunctions;
		}


		template<class Shape, class ShapeFunction>
		const int CFiniteElement<Shape, ShapeFunction, bool, bool>::GetNumberOfNodes() const
		{
			return m_shape.GetNumberOfNodes();
		}

		template<class Shape, class ShapeFunction>
		const int CFiniteElement<Shape, ShapeFunction, bool, bool>::GetDoFs() const
		{
			return m_shapefunctions.GetNumberOfShapeFunctions();
		}

		template<class Shape, class ShapeFunction>
		void CFiniteElement<Shape, ShapeFunction, bool, bool>::SetNeighbour(const int k, const int elem)
		{
			m_neighbours[k] = elem;
		}

		template<class Shape, class ShapeFunction>
		void CFiniteElement<Shape, ShapeFunction, bool, bool>::SetShapeFunction(const int k, const ShapeFunction& func)
		{
			m_shapefunctions = func;
		}

		template<class Shape, class ShapeFunction>
		void CFiniteElement<Shape, ShapeFunction, bool, bool>::SetShape(const Shape &shape)
		{
			m_shape = shape;
		}

		template<class Shape, class ShapeFunction>
		void CFiniteElement<Shape, ShapeFunction, bool, bool>::SetType(const int k)
		{
			m_type = k;
		}

		template<class Shape, class ShapeFunction>
		void CFiniteElement<Shape, ShapeFunction, bool, bool>::SetNode(const int k, const int node)
		{
			m_shape.SetNode(k, node);
		}
		template<class Shape, class ShapeFunction>
		const double CFiniteElement<Shape, ShapeFunction, bool, bool>::GetShapeFunction(const int k, const Mesh::Point &p) const
		{
			return m_shapefunctions.GetShapeFunction(k, p);
		}

		template<class Shape, class ShapeFunction>
		const Point CFiniteElement<Shape, ShapeFunction, bool, bool>::GetGradShapeFunction(const int k, const Mesh::Point &p) const
		{
			return m_shapefunctions.GetGradShapeFunction(k, p);
		}
		
		template<class Shape, class ShapeFunction>
		const int CFiniteElement<Shape, ShapeFunction, bool, bool>::IncreaseOrder()
		{
			if(m_shape.IncreaseOrder())
				return 1;
			if (m_shapefunctions.IncreaseOrder())
				return 1;
			return 0;
		}
		template<class Shape, class ShapeFunction>
		const double CFiniteElement<Shape, ShapeFunction, bool, bool>::GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
		{
			return m_shapefunctions.GetWeight(node, verts, f);
		}
		// fin.

	}

}


#endif /* CORENC_MESH_FINITEELEMENT_H_ */
