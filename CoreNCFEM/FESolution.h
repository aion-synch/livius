#ifndef CORENC_FESOLUTION_H
#define CORENC_FESOLUTION_H
/**
	Usually it is a vector<double> but some methods required different types like vector of vectors or double/tripple values
	The interface for dealing with solutions;
    NOT IN USE
**/

#include <vector>
#include "Point.h"
namespace corenc
{
	class CSolution
	{
	public:
		CSolution() {};
		virtual ~CSolution() {};
	};
	
	class CFESolution :public CSolution
	{
	public:
		CFESolution() :m_w{ 0 } {};
		~CFESolution() {}
		CFESolution& operator=(const CFESolution& fe)
		{
			m_w = fe.m_w;
			return *this;
		}
		CFESolution& operator=(const double fe)
		{
			m_w = fe;
			return *this;
		}
		CFESolution(const CFESolution& fe) :m_w{ fe.m_w } {}
		CFESolution(const double& fe) : m_w{ fe } {}
		operator double() const { return m_w; }
		/*double& operator=(const double fe)
		{
			m_w = fe;
			return m_w;
		}*/
		const bool operator==(const CFESolution& fe)
		{
			if (fe.m_w == m_w)
				return true;
			return false;
		}
		const bool operator!=(const CFESolution& fe)
		{
			if (fe.m_w != m_w)
				return true;
			return false;
		}
		CFESolution& operator+=(const CFESolution& fe)
		{
			m_w += fe.m_w;
			return *this;
		}
		CFESolution& operator-=(const CFESolution& fe)
		{
			m_w -= fe.m_w;
			return *this;
		}
		CFESolution& operator*=(const CFESolution& fe)
		{
			m_w *= fe.m_w;
			return *this;
		}
		CFESolution& operator/=(const CFESolution& fe)
		{
			m_w /= fe.m_w;
			return *this;
		}
		friend const double operator*(const CFESolution& lhs, const CFESolution& rhs)
		{
			return lhs.m_w * rhs.m_w;
		}
		friend const double operator*(const CFESolution& lhs, const double rhs)
		{
			return lhs.m_w * rhs;
		}
		friend const double operator*(const double lhs, const CFESolution& rhs)
		{
			return lhs * rhs.m_w;
		}
		friend const double operator-(const CFESolution& lhs, const CFESolution& rhs)
		{
			return lhs.m_w - rhs.m_w;
		}
		friend const double operator+(const CFESolution& lhs, const CFESolution& rhs)
		{
			return lhs.m_w + rhs.m_w;
		}
		friend const double operator/(const CFESolution& lhs, const CFESolution& rhs)
		{
			return lhs.m_w / rhs.m_w;
		}
	private:
		double m_w;
	};
	class CVecSolution :public CSolution
	{
	public:
		CVecSolution() :m_w{ 0 } {};
		~CVecSolution() {}
		std::vector<double> m_w;
	};
	
	class CFEweights
	{
	public:
		CFEweights() {};
		~CFEweights()
		{
			if (m_w.size() > 0)
				std::vector<CFESolution>().swap(m_w);
		};
		const CFESolution				getWeight(const unsigned int i) const { return m_w[i]; };
		const int						updateWeight(const unsigned int i, const CFESolution& cfe) 
		{ 
			if (i < m_w.size())
			{
				m_w[i] = cfe;
				return 0;
			}
			return 1;
		}
	private:
		std::vector<CFESolution>		m_w;
	};
}

#endif // !CORENC_FESOLUTION_H

