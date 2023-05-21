#pragma once
#ifndef CORENC_MESH_Point_h
#define CORENC_MESH_Point_h
#include <cmath>
#include <vector>
namespace corenc
{
	enum class Terms
	{
		// left-side
		// uv
		IUV,
		// grad u grad v
		IDUDV,
		// grad u v
		IDUV,
		// u grad v
		IUDV,
		// right-side
		EUV,
		EDUDV,
		EDUV,
		EUDV,
		EFV,
		// right-side matrix
		RUV,
        SUPG,
	};
	
	enum class Parameters
	{
		DIFFUSION,
		MASS,
		ADVECTION
	};
	
	namespace Mesh
	{
		class Point
		{
		public:
			Point() :x{ 0 }, y{ 0 }, z{ 0 } {}
			Point(const double _x, const double _y) :
				x{ _x }, y{ _y }, z{ 0 } {}
			Point(const double _x, const double _y, const double _z) :
				x{ _x }, y{ _y }, z{ _z } {}
			Point(const Point& p) :
				x{ p.x }, y{ p.y }, z{ p.z } {}
			double x, y, z;
			const double Jacobian() const { return 1; }
			Point& operator=(const Point& p)
			{
				x = p.x;
				y = p.y;
				z = p.z;
				return *this;
			}
			const bool operator==(const Point& p)
			{
				const double eps{ 1e-13 };
				if (fabs(x - p.x) < eps)
					if (fabs(y - p.y) < eps)
						if (fabs(z - p.z) < eps)
							return true;
				return false;
			}
			friend const bool operator!=(const Point& p1, const Point& p2)
			{
				const double eps{ 1e-13 };
				if (fabs(p1.x - p2.x) < eps)
					if (fabs(p1.y - p2.y) < eps)
						if (fabs(p1.z - p2.z) < eps)
							return false;
				return true;
			}
			const bool operator<(const Point& p2)
			{
				return (x < p2.x);
			}
			friend const double operator*(const Point& lhs, const Point& rhs)
			{
				return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
			}
			const Point operator*(const double rhs)
			{
				return Point{ x * rhs, y * rhs, z * rhs };
			}
			Point& operator+=(const Point& rhs)
			{
				x += rhs.x;
				y += rhs.y;
				z += rhs.z;
				return *this;
			}
			Point& operator*=(const double rhs)
			{
				x *= rhs;
				y *= rhs;
				z *= rhs;
				return *this;
			}
			friend const Point operator*(const Point& lhs, const double rhs)
			{
				return Point{ rhs * lhs.x, rhs * lhs.y, rhs * lhs.z };
			}
			friend const Point operator*(const double lhs, const Point& rhs)
			{
				return Point{ lhs * rhs.x, lhs * rhs.y, lhs * rhs.z };
			}
			friend const Point operator+(const Point& lhs, const Point& rhs)
			{
				return Point{ lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
			}
			friend const Point operator-(const Point& lhs, const Point& rhs)
			{
				return Point{ lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
			}
		};

		struct GaussTriangle
		{
			const static double m_tra[];
			const static double m_trb[];
			const static double m_sqrt15;
			const static double m_trw[];
			const static int	m_order;

		};
		struct GaussRectangular
		{
			const static double m_ra[];
			const static double m_rb[];
			const static double m_rw[];
			const static double m_a;
			const static double m_b;
			const static double m_c;
			const static double m_wa;
			const static double m_wb;
			const static double m_wc;
		};
		struct Gauss1dim
		{
			const static int	m_order;
			const static double m_a[];
			const static double m_sqrt35;
			const static double m_w[];
		};

		template<int N>
		struct Gauss1dimN
		{
			const static int	m_order;
			const static double m_a[];
			const static double m_w[];
		};


		struct GaussTetrahedron
		{
			const    static double    m_la[];
			const    static double    m_lb[];
			const    static double    m_lc[];
			const    static double    m_ld[];
			const    static double    m_w[];
			const    static double    m_psq, m_msq;
		};

		struct GaussRectangularCubic
		{
			const static double m_ra[];
			const static double m_rb[];
			const static double m_rc[];
			const static double m_rw[];
			const static double m_a;
			const static double m_b;
			const static double m_c;
			const static double m_w1;
			const static double m_w2;
			const static double m_w3;
			const static double m_w4;
			const static int m_s{ 34 };
		};
	}
}
#endif /* CORENC_MESH_Point_h */
