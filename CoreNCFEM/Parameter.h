// OK. DESCRIPTION.
// Here the known parameters are described. it is used then with meshes and problems etc.
#pragma once
#ifndef CORENC_MESH_PARAMETER_H_
#define CORENC_MESH_PARAMETER_H_

#include "Point.h"
#include <functional>
namespace corenc
{
	namespace Mesh
	{
		template<class T>
		class parameter
		{
		public:
			using cfunc = std::function<const T(const int, const int, const Point&)>;
			using cfunc_old = std::function<const T(const int, const Point&)>;
			parameter() :m_func{ [=](const int, const int, const Point&) {return T(); } } {};
			parameter(const cfunc& func):m_func{func}{}
			parameter(const cfunc_old& func)
			{ 
				cfunc f = [=](const int, const int n, const Point& p) {return func(n, p); };
				m_func = f;
			}
			parameter(const double _p) :m_func{ [=](const int, const int, const Point&) {return _p; } } {}
			parameter(const Mesh::Point _p) :m_func{ [=](const int, const int, const Point&) {return _p; } } {}
			parameter(const parameter<T>& _p) :m_func{ _p.m_func } {}
			~parameter() {};
			const T			get(const Point& p) const { return m_func(0, 0, p); };
			const T			get(const int number, const Point& p) const { return m_func(0, number, p); };
			const T			get(const int element, const int node, const Point& p) const { return m_func(element, node, p); };
			void			set(const cfunc& func) { m_func = func; };
		private:
			cfunc			m_func;
		};

        template<class T>
        class point_source
        {
        public:
            point_source() : m_point(Mesh::Point(0,0,0)), m_value(T(0)) {};
            point_source(const Mesh::Point& p, const T& val) : m_point(p), m_value(val) {};
            const T             get_value() const { return m_value; };
            const Mesh::Point   get_point() const { return m_point; };
            point_source<T>&    operator=(const point_source<T>& ps)
            {
                m_point = ps.m_point;
                m_value = ps.m_value;
                return *this;
            }
        private:
            Mesh::Point     m_point;
            T               m_value;
        };
		class CParameter
		{
		public:
			CParameter();
			//CParameter(const double _diff, const double _adv, const double _mass);
			CParameter(const parameter<double>& _diff, const parameter<double>& _adv, const parameter<double>& _mass);
			CParameter(const Parameters&, const parameter<double>&);
			~CParameter();
			const double		GetDiffusion() const;
			const double		GetAdvection() const;
			const double		GetMass() const;
			const double		GetDiffusion(const Point&) const;
			const double		GetAdvection(const Point&) const;
			const double		GetMass(const Point&) const;
		private:
			parameter<double>	m_diffusion;
			parameter<double>	m_advection;
			parameter<double>	m_mass;
		};
	}
}

#endif // !CORENC_MESH_PARAMETER_H_
