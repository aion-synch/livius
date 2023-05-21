#ifndef CORENC_METHODS_RUNGEKUTTA
#define CORENC_METHODS_RUNGEKUTTA

#include <memory>
#include "../Point.h"
#include <functional>
namespace corenc
{
	namespace method
	{
		template<class Problem, class Type>
		class RungeKutta
		{
		public:
			RungeKutta() {};
			RungeKutta(const double step, const double final, Problem* problem, const Type* solution) :
				m_step{ step }, 
				m_final{ final },
				m_problem{problem} {}
			const Type				discretize(const Type& solution, const std::function<const Type(const double time, const double time_step, const Type& curr_sol, Type* result)>& func);
			const Type				explicitEuler(const Type& solution, const std::function<const Type(const double time, const double time_step, const Type& curr_sol, Type* result)>& func);
			void					updateTimestep(const double step) { m_step = step; };
			~RungeKutta() {};
		private:
			double					m_step;
			double					m_final;
			double					m_curr;
			Problem*				m_problem;
			Type*					m_solution;
			static const std::vector<double> vector_mult(const std::vector<double>& lhs, const double rhs)
			{
				std::vector<double> vc(lhs);
				for (auto &it : vc)
					it *= rhs;
				return vc;
			}

			static const std::vector<double> vector_mult(const double rhs, const std::vector<double>& lhs)
			{
				std::vector<double> vc(lhs);
				for (auto &it : vc)
					it *= rhs;
				return vc;
			}

			static const std::vector<double> vector_divide(const std::vector<double>& lhs, const double rhs)
			{
				std::vector<double> vc(lhs);
				for (auto &it : vc)
					it /= rhs;
				return vc;
			}

			static const std::vector<double> vector_divide(const double rhs, const std::vector<double>& lhs)
			{
				std::vector<double> vc(lhs);
				for (auto &it : vc)
					it /= rhs;
				return vc;
			}

			static const std::vector<double> vector_add(const std::vector<double>& rhs, const std::vector<double>& lhs)
			{
				std::vector<double> vc(lhs);
				for (unsigned i{ 0 }; i < vc.size(); ++i)
					vc[i] += rhs[i];
				return vc;
			}
		};

		

		template<class Problem, class Type>
		const Type RungeKutta<Problem, Type>::discretize(const Type& u_pr, const std::function<const Type(const double time, const double time_step, const Type& curr_sol, Type* result)>& func)
		{
			Type k[4];
			const int n{ int(m_final / m_step) };
			func(m_curr, m_step, u_pr, &k[0]);
			//std::vector<double> tempc(m_curr.size());
			std::vector<double> tempu(u_pr.size());
			std::vector<double> tempk(u_pr.size());
			tempk = vector_divide(k[0], 2);
			tempu = vector_add(u_pr, tempk);
			func(m_curr + m_step / 2, m_step, tempu, &k[1]);
			//func(m_curr + m_step / 2, m_step, u_pr + k[0] / 2, &k[1]);

			tempk = vector_divide(k[1], 2);
			tempu = vector_add(u_pr, tempk);
			func(m_curr + m_step / 2, m_step, tempu, &k[2]);
			//func(m_curr + m_step / 2, m_step, u_pr + k[1] / 2, &k[2]);

			tempu = vector_add(u_pr, k[2]);
			func(m_curr + m_step, m_step, tempu, &k[3]);
			//func(m_curr + m_step, m_step, u_pr + k[2], &k[3]);

			tempk = vector_mult(k[1], 2);
			tempu = vector_mult(k[2], 2);
			k[3] = vector_add(k[3], tempu);
			k[3] = vector_add(k[3], tempk);
			k[3] = vector_add(k[3], k[0]);
			k[3] = vector_divide(k[3], 6.);
			//k[3] = k[0] + 2 * k[1] + 2 * k[2] + k[3];
			//k[3] = 1. / 6 * k[3];
			m_problem->addTerm(Terms::EUV);
			m_problem->addTerm(Terms::IUV);
			m_curr += m_step;
			return k[3];
		}
		
		template<class Problem, class Type>
		const Type RungeKutta<Problem, Type>::explicitEuler(const Type& u_pr, const std::function<const Type(const double time, const double time_step, const Type& curr_sol, Type* result)>& func)
		{
			Type k;
			func(m_curr, m_step, u_pr, &k);
			m_problem->addTerm(Terms::EUV);
			m_problem->addTerm(Terms::IUV);
			m_curr += m_step;
			return k;
		}
	}
}
#endif // !CORENC_METHODS_RUNGEKUTTA
