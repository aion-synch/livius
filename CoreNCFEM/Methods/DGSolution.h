#ifndef CORENC_METHODS_DGSOLUTION_H_
#define CORENC_METHODS_DGSOLUTION_H_

#include "DGMethod.h"
namespace corenc
{
	namespace method
	{
		template<class Grid>
		class DGSolution
		{
		public:
			DGSolution() {};
			DGSolution(const std::vector<double>& w) :m_w{ w } {}
			DGSolution(const DGSolution<Grid>& dg) :m_w{ dg.m_w } {}
			DGSolution<Grid>& operator=(const DGSolution<Grid>& dg)
			{
				m_w = dg.m_w;
				return *this;
			}
			~DGSolution()
			{
				if (m_w.size() > 0)
					std::vector<double>().swap(m_w);
			}
			const double				getWeight(const Grid& g, const Mesh::Point& p) const
			{
				if (m_w.size() > 0)
					return DGMethod<int, Grid, int>::GetSolution(g, m_w, p);
				return 0.;
			};
			const std::vector<double>	getWeights() const { return m_w; }
			const int					updateWeight(const unsigned int i, const double val)
			{
				if (i < m_w.size())
				{
					m_w[i] = val;
					return 0;
				}
				return 1;
			}
		private:
			std::vector<double>			m_w;
		};
		template<class Grid>
		class STSolution
		{
		public:
			STSolution() {};
			STSolution(const Grid& g):m_grid{g}{}
			STSolution(
				const std::vector<DGSolution<Grid>>& w,
				const std::vector<double> time,
				const Grid& g) : m_w{ w }, m_time{ time }, m_grid{g} {}
			STSolution(const STSolution<Grid>& st) :m_w{ st.m_w }, m_time{ st.m_time }, m_grid{st.m_grid} {}
			STSolution<Grid>& operator=(const STSolution<Grid>& st)
			{
				m_w = st.m_w;
				m_time = st.m_time;
				m_grid = st.m_grid;
				return *this;
			}
			~STSolution()
			{
				if (m_w.size() > 0)
					std::vector<DGSolution<Grid>>().swap(m_w);
				if (m_time.size() > 0)
					std::vector<double>().swap(m_time);
			}
			const double				getWeight(const Mesh::Point& p, const double time) const
			{
				int i = 0;
				auto sz = m_time.size();
				if (fabs(time) < 1e-14)
					return DGMethod<Grid>::GetSolution(m_grid, m_w[0].getWeights(), p);
				for (; i < sz; ++i)
				{
					if (time < m_time[i])
						break;
				}
				if (i == sz)
					--i;
				double dt = m_time[i] - m_time[i - 1];
				auto temp = DGMethod<Grid>::GetSolution(m_grid, m_w[i - 1].getWeights(), p);
				double du = DGMethod<Grid>::GetSolution(m_grid, m_w[i].getWeights(), p) - temp;
				return temp + du * (time - m_time[i - 1]) / dt;
			};
			const int					updateWeight(
				const std::vector<double> time,
				const std::vector<DGSolution<Grid>> w
			)
			{
				m_time = time;
				m_w = w;
			}
			const int					addTimeLayer(
				const double time,
				const DGSolution<Grid> w
			)
			{
				m_time.push_back(time);
				m_w.push_back(w);
				return 0;
			}
			const std::vector<DGSolution<Grid>>		getWeights() const { return m_w; }
		private:
			std::vector<DGSolution<Grid>>							m_w;
			std::vector<double>										m_time;
			Grid													m_grid;
		};
	}
}
#endif // !CORENC_METHODS_DGSOLUTION_H_
