// NO GENERALIZATION HERE
// JUST PLAIN DG FOR SYSTEM IN N - DIMENSIONAL SPACE FOR ONE TIME STEP
// CONSTANT BASIS FUNCTIONS

#pragma once
#ifndef CORENC_METHODS_SYSTEM_DG_METHOD_H_
#define CORENC_METHODS_SYSTEM_DG_METHOD_H_
#include <functional>
#include <set>
#include "../Point.h"
#include <memory>
#include <cmath>
#include "FEMethod.h"
#include <map>
#include <algorithm>
#include <vector>
#include "dg_flux.h"

namespace corenc
{
	namespace method
	{
		template<class Problem, class Grid, class Matrix>
		class system_dg_method
		{
		public:
			system_dg_method() :
				m_problem{ nullptr },
				m_CoarseGrid{ nullptr },
				m_GlobalMatrix{ nullptr },
				m_rhsvector{ nullptr }
			{};
			system_dg_method(
				Problem* p,
				Grid* g,
				Matrix* m,
				//Solution* s,
				const size_t sys_size,
				std::vector<double>* rhs):
				//const std::function<const double(const double)>& flux_function,
				//const DGFlux flux_type) :
				m_problem{ p },
				m_CoarseGrid{ g },
				m_GlobalMatrix{ m },
				m_N{ g->GetNumberOfElements() },
				m_Ns{ g->GetNumberOfBoundaries() },
				m_rhsvector{ rhs },
				//m_flux(flux_function),
				m_sys_size{sys_size}{
				GeneratePortrait();
			}
			~system_dg_method() {};
			const int					Assemble();
			const int					changeFlux(const DGFlux flux_type) { m_fluxtype = flux_type; return 0; };
			const Matrix*				GetGlobalMatrix() const { return m_GlobalMatrix; };
			const std::vector<double>	GetSolution() const { return m_vec; };
			const double				GetSolution(const std::vector<double>& point) const;
			const double				GetMaxSolution() const;
			const double				GetMinSolution() const;
			static const double	 		GetSolution(const Grid& g, const std::vector<double> &dg_sol, const Mesh::Point& p)
			{
				double sum{ 0 };
				auto nfem{ g.FindElement(p) };
				auto elem{ g.GetElement(nfem) };
				auto dofs{ elem->GetDoFs() };
				for (auto i{ 0 }; i < dofs; ++i)
				{
					sum += dg_sol[nfem * dofs + i] * elem->GetShapeFunction(i, p);
				}
				return sum;
			}
			const double	 		GetSolution(const std::vector<double> &dg_sol, const Mesh::Point& p)
			{
				double sum{ 0 };
				auto nfem{ m_CoarseGrid->FindElement(p) };
				auto elem{ m_CoarseGrid->GetElement(nfem) };
				auto dofs{ elem->GetDoFs() };
				for (auto i{ 0 }; i < dofs; ++i)
				{
					sum += dg_sol[nfem * dofs + i] * elem->GetShapeFunction(i, p);
				}
				return sum;
			}
			const int					toDGSolution(const Grid& g, std::vector<double>& dg_result) const
			{
				//dg_result->resize(m_rhsvector->size());
				dg_result.resize(m_rhsvector->size());
				for (unsigned i{ 0 }; i < g.GetNumberOfElements(); ++i)
				{
					auto elem{ g.GetElement(i) };
					auto dofs{ elem->GetDoFs() };
					for (unsigned j{ 0 }; j < dofs; ++j)
						//dg_result->operator[](m_nums[i] + j) = g.getSolution(i, j);
						dg_result[m_nums[i] + j] = g.getSolution(i, j);
				}
				return 0;
			}
			const int					updateWeights(const std::vector<double>& dg_result)
			{
				for (unsigned int i{ 0 }; i < (unsigned int)m_CoarseGrid->GetNumberOfElements(); ++i)
				{
					for (unsigned int j{ 0 }; j < (unsigned int)m_CoarseGrid->GetElement(i)->GetDoFs(); ++i)
						m_CoarseGrid->updateSolution(i, j, dg_result[m_nums[i] + j]);
				}
				return 0;
			}

			const int					DGtostandart(const std::vector<double>& dg_result)
			{
				for (unsigned int i{ 0 }; i < (unsigned int)m_CoarseGrid->GetNumberOfElements(); ++i)
				{
					auto elem{ m_CoarseGrid->GetElement(i) };
					auto dofs{ elem->GetDoFs() };
					for (unsigned int j{ 0 }; j < (unsigned int)dofs; ++j)
						//m_CoarseGrid->updateSolution(i, j, dg_result[m_nums[i] + j]);
						m_CoarseGrid->updateSolution(i, j, dg_result[m_nums[i] + j]);
				}
				return 0;
			}
		private:
			const int					GeneratePortrait();
			void						assembleBoundaries();
			void						assemble_flux(const unsigned boundary);
			const double				numerical_flux(const double ul, const double ur, const double fl, const double fr) const;
			void						MainConditions();
			void						SecondConditions();
			void						ThirdConditions();
			const int					AssembleGlobal();
			const int					AssembleFluxMatrix();
			Grid*						m_CoarseGrid;
			Matrix*						m_GlobalMatrix;
			std::vector<double>*		m_rhsvector;
			std::vector<unsigned int>	m_nums;
			unsigned int				m_N;  // number of elements
			unsigned int				m_Ns; // number of boundaries
			unsigned int				m_size;
			Problem*					m_problem;
			std::vector<double>			m_vec;
			std::vector<double>			m_solution;
			//std::function<const Mesh::Point(const Mesh::Point)> m_numflux;
			//std::function<const Mesh::Point(const Mesh::Point)> m_flux;
			DGFlux						m_fluxtype;
			size_t						m_sys_size;
			std::function<const double(const double)> m_flux;
			const int					AssembleLocalMatrix(const int);
		};

		template<class Grid>
		class system_dg_method<Grid, bool, bool>
		{
		public:
			static const double	 		GetSolution(const Grid& g, const std::vector<double> &dg_sol, const Mesh::Point& p)
			{
				double sum{ 0 };
				auto nfem{ g.FindElement(p) };
				auto elem{ g.GetElement(nfem) };
				auto dofs{ elem->GetDoFs() };
				for (auto i{ 0 }; i < dofs; ++i)
				{
					sum += dg_sol[nfem * dofs + i] * elem->GetShapeFunction(i, p);
				}
				return sum;
			}
		};

		template<class Problem, class  Grid, class Matrix>
		const int system_dg_method<Problem, Grid, Matrix>::Assemble()
		{
			//GeneratePortrait();
			AssembleGlobal();
			AssembleFluxMatrix();
			MainConditions();
			SecondConditions();
			ThirdConditions();
			return 0;
		}
		template<class Problem, class Grid, class Matrix>
		const int system_dg_method<Problem, Grid, Matrix>::GeneratePortrait()
		{
			int lorder, rorder, order;
			std::vector<std::set<unsigned int>> temp;
			unsigned int i, j, nk, ne, k, sz, size;
			m_size = 0;
			m_nums.resize(m_N * m_sys_size);

			nk = 0;
			sz = m_N * m_sys_size;
			for (i = 0, k = 0; k < sz; ++i, k += m_sys_size)
			{
				size = m_CoarseGrid->GetElement(i)->GetDoFs() * m_sys_size;
				for(j = 0; j < m_sys_size; ++j)
					m_nums[k + j] = m_size + j * m_CoarseGrid->GetElement(i)->GetDoFs();
				m_size += size;
			}
			temp.resize(m_size);
			sz = m_Ns;
			for (k = 0; k < sz; k += m_sys_size)
			{
				auto bound = m_CoarseGrid->GetBoundary(k);
				nk = bound->GetNeighbour(0);
				ne = bound->GetNeighbour(1);
				lorder = m_CoarseGrid->GetElement(nk)->GetDoFs();
				if (ne != -1)
				{

					rorder = m_CoarseGrid->GetElement(ne)->GetDoFs();
					for (i = 0; i < lorder; ++i)
						for (j = 0; j < rorder; ++j)
							temp[m_nums[ne] + j].insert(m_nums[nk] + i);
					for (i = 0; i < lorder; ++i)
						for (j = i + 1; j < lorder; ++j)
							temp[m_nums[nk] + j].insert(m_nums[nk] + i);
				}
				else
				{
					for (i = 0; i < lorder; ++i)
						for (j = i + 1; j < lorder; ++j)
							temp[m_nums[nk] + j].insert(m_nums[nk] + i);
				}
			}

			/*temp.resize(m_CoarseGrid->GetNumberOfNodes());
			m_nums.resize(m_CoarseGrid->GetNumberOfNodes());
			lorder = m_CoarseGrid->GetElement(0)->GetDoFs();
			for (k = 0; k < m_CoarseGrid->GetNumberOfNodes(); ++k)
			m_nums[k] = k;
			//for (auto elem : m_CoarseGrid->GetElements())
			for(k = 0; k < m_CoarseGrid->GetNumberOfElements(); ++k)
			{
			auto elem{ m_CoarseGrid->GetElement(k) };
			auto order{ elem->GetDoFs() };
			for (i = 0; i < order; ++i)
			for (j = 0; j < order; ++j)
			if (elem->GetNode(j) > elem->GetNode(i))
			temp[elem->GetNode(j)].insert(elem->GetNode(i));
			}*/
			m_GlobalMatrix->Create(temp.size(), temp);
			m_rhsvector->resize(temp.size());
			//m_vec.resize(temp.size());
			return 0;
		}

		template<class Problem, class Grid, class Matrix>
		const int system_dg_method<Problem, Grid, Matrix>::AssembleLocalMatrix(const int l)
		{
			int i, j, k, nodes;
			double mij;
			const auto& elem{ m_CoarseGrid->GetElement(l) };
			const auto& dofs{ elem->GetDoFs() };
			nodes = elem->GetNumberOfNodes();
			std::vector<Mesh::Point> points(nodes);
			for (i = 0; i < nodes; ++i)
				points[i] = m_CoarseGrid->GetNode(elem->GetNode(i));
			for (k = 0; k < m_problem->getNumberOfTerms(); ++k)
			{
				switch (m_problem->getTerm(k))
				{
				case Terms::EUV:
					//for (i = 0; i < dofs; ++i)
					//{
					//	for (j = 0; j < dofs; ++j)
					//	{
					//		auto M = [&](const Mesh::Point& p)
					//		{
					//			return elem->GetShapeFunction(i, p) * elem->GetShapeFunction(j, p);
					//		};
					//		mij = elem->Integrate(M, points);
					//		m_rhsvector->operator[](m_nums[l] + i) += m_CoarseGrid->getParameter(Parameters::MASS, l, j) * m_CoarseGrid->getSolution(l, j) * mij;
					//	}
					//}
					for(size_t j = 0; j < m_sys_size; ++j)
						m_rhsvector->operator[](m_nums[l] + j) += m_problem->get_solution(j, l, elem->GetType(), points[l]);
					break;
				default:
					break;
				}
			}
			return 0;
		}

		template<class Problem, class Grid, class Matrix>
		const int system_dg_method<Problem, Grid, Matrix>::AssembleGlobal()
		{
			for (int l = 0; l < m_N; ++l)
				AssembleLocalMatrix(l);
			return 0;
		}

		template<class Problem, class Grid, class Matrix>
		const int system_dg_method<Problem, Grid, Matrix>::AssembleFluxMatrix()
		{
			auto Nb{ m_CoarseGrid->GetNumberOfBoundaries() };
			unsigned int l;
			switch (m_fluxtype)
			{
			case corenc::method::DGFlux::ELaxFriedrichs:

				for (l = 0; l < Nb; ++l)
				{
					const auto& bound{ m_CoarseGrid->GetBoundary(l) };
					const auto& nk{ bound->GetNeighbour(0) };
					const auto& ne{ bound->GetNeighbour(1) };
					const auto& elemk{ m_CoarseGrid->GetElement(nk) };
					const auto& dofs{ bound->GetDoFs() };
					const auto& dofsk{ elemk->GetDoFs() };
					double C{ 0 };
					unsigned int i, j;
					std::vector<Mesh::Point> points(dofs);
					for (i = 0; i < dofs; ++i)
						points[i] = m_CoarseGrid->GetNode(bound->GetNode(i));
					if (ne > -1)
					{
						const auto& eleme{ m_CoarseGrid->GetElement(ne) };
						const auto& dofse{ eleme->GetDoFs() };
						for (i = 0; i < dofsk; ++i)
						{
							for (j = 0; j < dofsk; ++j)
							{
								auto Mkk = [&](const Mesh::Point& p)
								{
									return elemk->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
								};
								auto temp{ bound->Integrate(Mkk, points) };
								C = std::max(fabs(m_CoarseGrid->getSolution(ne, i)), fabs(m_CoarseGrid->getSolution(nk, j)));
								//m_rhsvector->operator[](m_nums[nk] + i) += -0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) * temp - C * m_CoarseGrid->getSolution(nk, j) * temp);
								auto val{ -0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) + C * m_CoarseGrid->getSolution(nk, j)) * temp };
								m_rhsvector->operator[](m_nums[nk] + i) += val;
								///kk[m_nums[nk] + i] += val;
								//lv[m_nums[nk] + i] += val;
							}
						}
						for (i = 0; i < dofsk; ++i)
						{
							for (j = 0; j < dofse; ++j)
							{
								auto Mke = [&](const Mesh::Point& p)
								{
									return eleme->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
								};
								auto temp{ bound->Integrate(Mke, points) };
								C = std::max(fabs(m_CoarseGrid->getSolution(nk, i)), fabs(m_CoarseGrid->getSolution(ne, j)));
								//m_rhsvector->operator[](m_nums[nk] +7 i) += -0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) * temp - C * m_CoarseGrid->getSolution(ne, j) * temp);
								auto val{ -0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) - C * m_CoarseGrid->getSolution(ne, j)) * temp };
								m_rhsvector->operator[](m_nums[nk] + i) += val;
								//ke[m_nums[nk] + i] += val;
								//lv[m_nums[nk] + i] += val;
							}
						}
						for (i = 0; i < dofse; ++i)
						{
							for (j = 0; j < dofsk; ++j)
							{
								auto Mek = [&](const Mesh::Point& p)
								{
									return eleme->GetShapeFunction(i, p) * elemk->GetShapeFunction(j, p);
								};
								auto temp{ bound->Integrate(Mek, points) };
								C = std::max(fabs(m_CoarseGrid->getSolution(nk, j)), fabs(m_CoarseGrid->getSolution(ne, i)));
								//m_rhsvector->operator[](m_nums[ne] + i) += 0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) * temp - C * m_CoarseGrid->getSolution(nk, j) * temp);
								auto val{ 0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) + C * m_CoarseGrid->getSolution(nk, j)) * temp };
								m_rhsvector->operator[](m_nums[ne] + i) += val;
								//ek[m_nums[ne] + i] += val;
								//rv[m_nums[ne] + i] += val;
							}
						}
						for (i = 0; i < dofse; ++i)
						{
							for (j = 0; j < dofse; ++j)
							{
								auto Mee = [&](const Mesh::Point& p)
								{
									return eleme->GetShapeFunction(j, p) * eleme->GetShapeFunction(i, p);
								};
								auto temp{ bound->Integrate(Mee, points) };
								C = std::max(fabs(m_CoarseGrid->getSolution(nk, i)), fabs(m_CoarseGrid->getSolution(ne, j)));
								//m_rhsvector->operator[](m_nums[ne] + i) += 0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) * temp - C * m_CoarseGrid->getSolution(ne, j) * temp);
								auto val{ 0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) - C * m_CoarseGrid->getSolution(ne, j)) * temp };
								m_rhsvector->operator[](m_nums[ne] + i) += val;
								//ee[m_nums[ne] + i] += val;
								//rv[m_nums[ne] + i] += val;
							}
						}
					}
					else
					{
						//C = m_flux(m_CoarseGrid->getSolution(nk, 0));
						//m_rhsvector->operator[](m_nums[nk]) = C;
						if (l == 0)
						{
							for (i = 0; i < dofsk; ++i)
							{
								for (j = 0; j < dofsk; ++j)
								{
									auto Mkk = [&](const Mesh::Point& p)
									{
										return elemk->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
									};
									auto temp{ bound->Integrate(Mkk, points) };
									//m_rhsvector->operator[](m_nums[nk]+i) -= ((ne+int(l))>0?-1:1)*(m_flux(m_CoarseGrid->getSolution(nk, j)) * temp - C * m_CoarseGrid->getSolution(nk, j) * temp);
									auto fl = m_flux(m_CoarseGrid->getSolution(nk, j)) * temp;
									m_rhsvector->operator[](m_nums[nk] + i) += fl * temp;
									//if(C >= 0)
									//m_rhsvector->operator[](m_nums[nk] + i) += ((ne + int(l))>0 ? 1 : 0) * C * temp;
									//m_rhsvector->operator[](m_nums[nk] + i) += 1e10 * C * temp;
								}
							}
						}
						else
						{
							for (i = 0; i < dofsk; ++i)
							{
								for (j = 0; j < dofsk; ++j)
								{
									auto Mkk = [&](const Mesh::Point& p)
									{
										return elemk->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
									};
									auto temp{ bound->Integrate(Mkk, points) };
									auto fl = m_flux(m_CoarseGrid->getSolution(nk, j)) * temp;
									m_rhsvector->operator[](m_nums[nk] + i) -= fl * temp;
								}
							}
						}
					}
				}

				// explicit LF flux
				break;
			default:
				break;
			}
			return 0;
		}

		template<class Problem, class Grid, class Matrix>
		void system_dg_method<Problem, Grid, Matrix>::assemble_flux(const unsigned l)
		{
			const auto& bound{ m_CoarseGrid->GetBoundary(l) };
			const auto& nk{ bound->GetNeighbour(0) };
			const auto& ne{ bound->GetNeighbour(1) };
			const auto& elemk{ m_CoarseGrid->GetElement(nk) };
			const auto& dofs{ bound->GetDoFs() };
			const auto& dofsk{ elemk->GetDoFs() };
			double C{ 0 };
			unsigned int i, j;
			std::vector<Mesh::Point> points(dofs);
			for (i = 0; i < dofs; ++i)
				points[i] = m_CoarseGrid->GetNode(bound->GetNode(i));
			C = 2;
			if (ne > -1)
			{
				const auto& eleme{ m_CoarseGrid->GetElement(ne) };
				const auto& dofse{ eleme->GetDoFs() };
				for (i = 0; i < dofsk; ++i)
				{
					for (j = 0; j < dofsk; ++j)
					{
						auto Mkk = [&](const Mesh::Point& p)
						{
							return elemk->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
						};
						auto temp{ bound->Integrate(Mkk, points) };
						C = std::max(m_CoarseGrid->getSolution(nk, i), m_CoarseGrid->getSolution(nk, j));
						auto val{ -0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) + C * m_CoarseGrid->getSolution(nk, j)) * temp };
						m_rhsvector->operator[](m_nums[nk] + i) += val;
						//kk[m_nums[nk] + i] += val;
						//lv[m_nums[nk] + i] += val;
					}
				}
				for (i = 0; i < dofsk; ++i)
				{
					for (j = 0; j < dofse; ++j)
					{
						auto Mke = [&](const Mesh::Point& p)
						{
							return eleme->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
						};
						auto temp{ bound->Integrate(Mke, points) };
						C = std::max(m_CoarseGrid->getSolution(nk, i), m_CoarseGrid->getSolution(ne, j));
						//m_rhsvector->operator[](m_nums[nk] +7 i) += -0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) * temp - C * m_CoarseGrid->getSolution(ne, j) * temp);
						auto val{ -0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) - C * m_CoarseGrid->getSolution(ne, j)) * temp };
						m_rhsvector->operator[](m_nums[nk] + i) += val;
						//ke[m_nums[nk] + i] += val;
						//rv[m_nums[nk] + i] += val;
					}
				}
				for (i = 0; i < dofse; ++i)
				{
					for (j = 0; j < dofsk; ++j)
					{
						auto Mek = [&](const Mesh::Point& p)
						{
							return eleme->GetShapeFunction(i, p) * elemk->GetShapeFunction(j, p);
						};
						auto temp{ bound->Integrate(Mek, points) };
						C = std::max(m_CoarseGrid->getSolution(nk, j), m_CoarseGrid->getSolution(ne, i));
						//m_rhsvector->operator[](m_nums[ne] + i) += 0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) * temp - C * m_CoarseGrid->getSolution(nk, j) * temp);
						auto val{ 0.5*(m_flux(m_CoarseGrid->getSolution(nk, j)) + C * m_CoarseGrid->getSolution(nk, j)) * temp };
						m_rhsvector->operator[](m_nums[ne] + i) += val;
						//ek[m_nums[ne] + i] += val;
						//lv[m_nums[ne] + i] += val;
					}
				}
				for (i = 0; i < dofse; ++i)
				{
					for (j = 0; j < dofse; ++j)
					{
						auto Mee = [&](const Mesh::Point& p)
						{
							return eleme->GetShapeFunction(j, p) * eleme->GetShapeFunction(i, p);
						};
						auto temp{ bound->Integrate(Mee, points) };
						C = std::max(m_CoarseGrid->getSolution(ne, i), m_CoarseGrid->getSolution(ne, j));
						//m_rhsvector->operator[](m_nums[ne] + i) += 0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) * temp - C * m_CoarseGrid->getSolution(ne, j) * temp);
						auto val{ 0.5*(m_flux(m_CoarseGrid->getSolution(ne, j)) - C * m_CoarseGrid->getSolution(ne, j)) * temp };
						m_rhsvector->operator[](m_nums[ne] + i) += val;
						//ee[m_nums[ne] + i] += val;
						//rv[m_nums[ne] + i] += val;
					}
				}
			}
			else
			{
				for (i = 0; i < dofsk; ++i)
				{
					for (j = 0; j < dofsk; ++j)
					{
						auto Mkk = [&](const Mesh::Point& p)
						{
							return elemk->GetShapeFunction(j, p) * elemk->GetShapeFunction(i, p);
						};
						auto temp{ bound->Integrate(Mkk, points) };
						auto fl = m_flux(m_CoarseGrid->getSolution(nk, j));
						C = m_CoarseGrid->getSolution(nk, j);
						m_rhsvector->operator[](m_nums[nk] + i) += ((ne + int(l))>0 ? 0 : 1) *fl * temp;
					}
				}
			}
		}
		template<class Problem, class Grid, class Matrix>
		void system_dg_method<Problem, Grid, Matrix>::MainConditions()
		{

		}

		template<class Problem, class Grid, class Matrix>
		void system_dg_method<Problem, Grid, Matrix>::SecondConditions()
		{

		}

		template<class Problem, class Grid, class Matrix>
		void system_dg_method<Problem, Grid, Matrix>::ThirdConditions()
		{

		}

		template<class Problem, class Grid, class Matrix>
		const double system_dg_method<Problem, Grid, Matrix>::GetMaxSolution() const
		{
			return 0.;
		}

		template<class Problem, class Grid, class Matrix>
		const double system_dg_method<Problem, Grid, Matrix>::GetMinSolution() const
		{
			return 0.;
		}

		template<class Problem, class Grid, class Matrix>
		const double system_dg_method<Problem, Grid, Matrix>::GetSolution(const std::vector<double>& point) const
		{
			return 0.;
		}
	}
}
#endif // !CORENC_METHODS_SYSTEM_DG_METHOD_H_
