#ifndef FINITEELEMENT2D_H
#define FINITEELEMENT2D_H
#include <functional>
#include <iostream>
#include <vector>
#include "../Point.h"
#include "../FESolution.h"
namespace corenc
{
    namespace Mesh
    {
    // 2d element

        template<class T = bool>
        class CElement2D;
        using function_dp = std::function<const double(const Point&)>;
        template<>
        class CElement2D<bool>
        {
        public:
            CElement2D() {}
            virtual ~CElement2D() {}
            virtual const int									GetType() const = 0;
            virtual CElement2D<>*								Clone() const = 0;
            virtual const int									GetDoFs() const = 0;
            virtual const int									GetNode(const int) const = 0;
            virtual const int									GetNeighbour(const int) const = 0;
            virtual void										SetNeighbour(const int k, const int elem) = 0;
            virtual void										SetType(const int) = 0;
            virtual void										SetNode(const int, const int) = 0;
            virtual const int									GetNumberOfNodes() const = 0;
            virtual const double								GetShapeFunction(const int, const Point&) const = 0;
            virtual const Point									GetGradShapeFunction(const int, const Point&) const = 0;
            virtual const Point									GetNormal() const = 0;
            virtual void										ReverseNormal() = 0;
            virtual const int                                   SetOrder(const int px, const int py) = 0;
            virtual const double								Integrate(const function_dp&, const std::vector<Point>& v) const = 0;
            virtual const Point									Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const = 0;
            virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const = 0;
            virtual const double								GetWeight(const int, const std::vector<Point>& verts, const function_dp& f) const = 0;
            virtual const int									IncreaseOrder() = 0;
            virtual const double								GetMeasure() const = 0;
        };

        template<class T>
        class CElement2D
        {
        public:
            CElement2D() {}
            virtual ~CElement2D() {}
            virtual const int									GetType() const = 0;
            virtual CElement2D*                                 Clone() const = 0;
            virtual const int									GetDoFs() const = 0;
            virtual const int									GetNode(const int) const = 0;
            virtual const int									GetNeighbour(const int) const = 0;
            virtual void										SetNeighbour(const int k, const int elem) = 0;
            virtual void										SetType(const int) = 0;
            virtual void										SetNode(const int, const int) = 0;
            virtual const int									GetNumberOfNodes() const = 0;
            virtual const double								GetShapeFunction(const int, const Point&) const = 0;
            virtual const Point									GetGradShapeFunction(const int, const Point&) const = 0;
            virtual const Point									GetNormal() const = 0;
            virtual void										ReverseNormal() = 0;
            virtual const int									IncreaseOrder() = 0;
            virtual const int                                   SetOrder(const int px, const int py) = 0;
            virtual const double								Integrate(const std::function<const double(const Point&)>&, const std::vector<Point>& v) const = 0;
            virtual const Point									Integrate(const std::function<const Point(const Point&)>&, const std::vector<Point>& v) const = 0;
            virtual const std::vector<double>					Integrate(const std::function<const std::vector<double>(const Point&)>&, const std::vector<Point>&) const = 0;
            virtual const double								GetMeasure() const = 0;
            virtual const double								GetWeight(const int, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const = 0;
        };

        template<class Shape, class ShapeFunction>
        class CFiniteElement2D : public CElement2D<>
        {
        public:
            CFiniteElement2D() {}
            CFiniteElement2D(const int* nodes, const Point* points, const int dofs) :
                m_shape{ nodes },
                m_shapefunctions{ points, dofs },
                m_type{ -1 } {
                m_neighbours[0] = -1; m_neighbours[1] = -1;
            };
            CFiniteElement2D(const int* nodes, const Point* points, const int dofs, const int type) :
                m_shape{ nodes, dofs },
                m_shapefunctions{ points, dofs },
                m_type{ type } {
                m_neighbours[0] = -1; m_neighbours[1] = -1;
            };
            CFiniteElement2D(const int* nodes, const Point* points) :
                m_shape{ nodes },
                m_shapefunctions{ points },
                m_type{ -1 } {
                m_neighbours[0] = -1; m_neighbours[1] = -1;
            };
            CFiniteElement2D(const Shape& shape, const ShapeFunction& f) :
                m_shape{ shape },
                m_shapefunctions{ f },
                m_type{ -1 } {
                m_neighbours[0] = -1; m_neighbours[1] = -1;
            };
            CFiniteElement2D(const Shape& shape, const ShapeFunction& shfunc, const int type) :
                m_shape{ shape },
                m_shapefunctions{ shfunc },
                m_type{ type } {
                m_neighbours[0] = -1; m_neighbours[1] = -1;
            };
            CFiniteElement2D(const Shape& shape, const ShapeFunction& shfunc, const int type, const int* neigs) :
                                m_shape{ shape },
                                m_shapefunctions{ shfunc },
                                m_type{ type } {
                                m_neighbours[0] = neigs[0]; m_neighbours[1] = neigs[1];
                        };
            CFiniteElement2D(const CFiniteElement2D&e) :
                m_shape{ e.m_shape },
                m_shapefunctions{ e.m_shapefunctions },
                m_type{ e.m_type } {
                m_neighbours[0] = e.m_neighbours[0]; m_neighbours[1] = e.m_neighbours[1];
            };
            friend const bool		operator==(const CFiniteElement2D& e1, const CFiniteElement2D& e2)
            {
                if (e1.m_shape == e2.m_shape)
                    return true;
                return false;
            }
            // don't forget to delete after the call
            CElement2D<>*			Clone() const
            {
                                return new CFiniteElement2D(m_shape, m_shapefunctions, m_type, m_neighbours);
            };
            ~CFiniteElement2D() {}
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
            const int               SetOrder(const int px, const int py);
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
            CFiniteElement2D&         operator=(const CFiniteElement2D& e)
            {
                m_shape = e.m_shape;
                m_shapefunctions = e.m_shapefunctions;
                m_type = e.m_type;
                return *this;
            }
            friend std::istream&	operator>>(std::istream& is, CFiniteElement2D& k)
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


        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::GetType() const
        {
            return m_type;
        }
        template<class Shape, class ShapeFunction>
        const double CFiniteElement2D<Shape, ShapeFunction>::GetMeasure() const
        {
            return m_shapefunctions.GetMeasure();
        }
        //template<class Shape, class ShapeFunction>
        //const int CFiniteElement2D<Shape, ShapeFunction>::SetValue(const int number, CSolution* value)
        //{
        //	return m_shapefunctions.SetValue(number, value);
        //}
        //template<class Shape, class ShapeFunction>
        //CSolution* CFiniteElement2D<Shape, ShapeFunction>::GetValue(const int p)
        //{
        //	auto&& val = m_shapefunctions.GetValue(p);
        //	return const_cast<CSolution*>(static_cast<const CSolution*>(&val));;
        //}
        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::GetNode(const int k) const
        {
            return m_shape.GetNode(k);
        }
        template<class Shape, class ShapeFunction>
        const double CFiniteElement2D<Shape, ShapeFunction>::Integrate(const std::function<const double(const Point&)>&f, const std::vector<Point>& v) const
        {
            return m_shape.Integrate(f, v);
        }

        template<class Shape, class ShapeFunction>
        const Point CFiniteElement2D<Shape, ShapeFunction>::Integrate(const std::function<const Point(const Point&)>&f, const std::vector<Point>& v) const
        {
            return m_shape.Integrate(f, v);
        }

        template<class Shape, class ShapeFunction>
        const std::vector<double> CFiniteElement2D<Shape, ShapeFunction>::Integrate(const std::function<const std::vector<double>(const Point&)>&f, const std::vector<Point>& v) const
        {
            return m_shape.Integrate(f, v);
        }
        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::GetNeighbour(const int k) const
        {
            return m_neighbours[k];
        }

        template<class Shape, class ShapeFunction>
        const Point CFiniteElement2D<Shape, ShapeFunction>::GetNormal() const
        {
            return m_shapefunctions.GetNormal();
        }

        template<class Shape, class ShapeFunction>
        void CFiniteElement2D<Shape, ShapeFunction>::ReverseNormal()
        {
            m_shapefunctions.ReverseNormal();
        }


        template<class Shape, class ShapeFunction>
        const Shape CFiniteElement2D<Shape, ShapeFunction>::GetShape() const
        {
            return m_shape;
        }

        template<class Shape, class ShapeFunction>
        const ShapeFunction CFiniteElement2D<Shape, ShapeFunction>::GetShapeFunctions() const
        {
            return m_shapefunctions;
        }


        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::GetNumberOfNodes() const
        {
            return m_shape.GetNumberOfNodes();
        }

        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::GetDoFs() const
        {
            return m_shapefunctions.GetNumberOfShapeFunctions();
        }

        template<class Shape, class ShapeFunction>
        void CFiniteElement2D<Shape, ShapeFunction>::SetNeighbour(const int k, const int elem)
        {
            m_neighbours[k] = elem;
        }

        template<class Shape, class ShapeFunction>
        void CFiniteElement2D<Shape, ShapeFunction>::SetShapeFunction(const int k, const ShapeFunction& func)
        {
            m_shapefunctions = func;
        }

        template<class Shape, class ShapeFunction>
        void CFiniteElement2D<Shape, ShapeFunction>::SetShape(const Shape &shape)
        {
            m_shape = shape;
        }

        template<class Shape, class ShapeFunction>
        void CFiniteElement2D<Shape, ShapeFunction>::SetType(const int k)
        {
            m_type = k;
        }

        template<class Shape, class ShapeFunction>
        void CFiniteElement2D<Shape, ShapeFunction>::SetNode(const int k, const int node)
        {
            m_shape.SetNode(k, node);
        }
        template<class Shape, class ShapeFunction>
        const double CFiniteElement2D<Shape, ShapeFunction>::GetShapeFunction(const int k, const Mesh::Point &p) const
        {
            return m_shapefunctions.GetShapeFunction(k, p);
        }

        template<class Shape, class ShapeFunction>
        const Point CFiniteElement2D<Shape, ShapeFunction>::GetGradShapeFunction(const int k, const Mesh::Point &p) const
        {
            return m_shapefunctions.GetGradShapeFunction(k, p);
        }

        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::IncreaseOrder()
        {
            if(m_shape.IncreaseOrder())
                return 1;
            if (m_shapefunctions.IncreaseOrder())
                return 1;
            return 0;
        }

        template<class Shape, class ShapeFunction>
        const int CFiniteElement2D<Shape, ShapeFunction>::SetOrder(const int px, const int py)
        {
            if(m_shape.SetOrder(px, py))
                return 1;
            if (m_shapefunctions.SetOrder(px, py))
                return 1;
            return 0;
        }

        template<class Shape, class ShapeFunction>
        const double CFiniteElement2D<Shape, ShapeFunction>::GetWeight(const int node, const std::vector<Point>& verts, const std::function<const double(const Point&)>& f) const
        {
            return m_shapefunctions.GetWeight(node, verts, f);
        }
        // fin.
    }
}

#endif // FINITEELEMENT2D_H
