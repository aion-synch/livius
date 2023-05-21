// FEMethod.h describes an abstract interface and functions for a general finite element method

#ifndef CORENC_METHODS_FEMethod_h
#define CORENC_METHODS_FEMethod_h
#include <functional>
#include <set>
#include "../Point.h"
#include "../Parameter.h"
#include "CSMethod.h"
#include <memory>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
namespace corenc
{
	namespace Mesh
	{
		class Point;
	}
	namespace method
	{
		enum class BoundaryType
		{
			MAIN,
			SECOND,
			THIRD,
			FREE
		};
		// class Type = Type of the solution, for ex vector or double, or even more specific

		
		template<class Type>
		class CFEMethod
		{
		public:
			CFEMethod() {};
			virtual ~CFEMethod() {};
			virtual const int							Assemble() = 0;
			virtual const Type							GetSolution(const std::vector<double>& point) const = 0;
			virtual const std::vector<Type>				GetSolution() const = 0;
			virtual const Type							GetMaxSolution() const = 0;
			virtual const Type							GetMinSolution() const = 0;
		};

		template<class Problem, class Grid, class Matrix>
		class FEMethod
		{
		public:
			FEMethod() :
				m_problem{nullptr},
				m_Grid{nullptr},
				m_GlobalMatrix{nullptr},
				m_RightMatrix{nullptr},
				m_rhsvector{nullptr}
			{}
			FEMethod(
				Problem* p,
				Grid* g,
				Matrix* m,
				std::vector<double>* rhs):
				m_problem{ p },
				m_Grid{ g->Clone() },
				m_GlobalMatrix{ m },
				m_N{ g->GetNumberOfElements() },
				m_Ns{ g->GetNumberOfBoundaries() },
				m_rhsvector{ rhs }{
				//GeneratePortrait();
			}
			FEMethod(
				Problem* p,
				Grid* g,
				Matrix* m,
				Matrix* rm,
				std::vector<double>* rhs):
				m_problem{ p },
				m_Grid{ g->Clone() },
				m_GlobalMatrix{ m },
				m_RightMatrix{ rm },
				m_N{ g->GetNumberOfElements() },
				m_Ns{ g->GetNumberOfBoundaries() },
				m_rhsvector{ rhs }{
				//GeneratePortrait();
			}
			FEMethod(const std::shared_ptr<Grid>& grid) :m_Grid{ grid->Clone() } {}
			FEMethod(Grid* grid) :m_Grid{ grid->Clone() } {}
			FEMethod(const FEMethod& meth) :
				m_Grid{ meth.m_Grid->Clone() },
				//m_GlobalMatrix{ meth.m_GlobalMatrix->Clone() },
				//m_rhsvector{ meth.m_rhsvector },
				//m_problem{ meth.m_problem },
				m_time{ meth.m_time },
				//m_solution{ meth.m_solution },
				m_size{ meth.m_size },
				m_N{ meth.m_N },
				m_Ns{ meth.m_Ns },
				m_nums{ meth.m_nums }
			{};
            FEMethod&                   operator=(const FEMethod& fem)
            {
                m_Grid = fem.m_Grid->Clone();
                m_time = fem.m_time;
                m_size = fem.m_size;
                m_N = fem.m_N;
                m_Ns = fem.m_Ns;
                m_nums = fem.m_nums;
                return *this;
            }
			void						Discretization();
			const double				GetValue(const Mesh::Point&) const;
			const double				GetValue(const Mesh::Point&, const std::vector<double>& vec) const;
			const double				GetValue(const Mesh::Point&, const std::vector<double>& vec, const int num) const;
			//const Mesh::Point			GetGradValue(const Mesh::Point&, const std::vector<double>& vec) const;
			//const Mesh::Point			GetLambdaGrad(const Mesh::Point&, const std::vector<double>& vec) const;
			const double				GetEffective(const std::vector<double>& vec) const;
			void						ProjectSolution(std::vector<double>&, std::function<const double(const Mesh::Point&, const std::vector<double>&, const int)> GetValue, std::vector<double>& sol);
			void						ProjectSolution(std::vector<double>&, std::function<const double(const Mesh::Point&, const std::vector<double>&)> GetValue, std::vector<double>& sol, const int);
			void						LoadSolution(const std::vector<double>& vec);
			const std::vector<double>	SetSolution(const int sol, const int liq, const double, const double, const double);
			void					 	GetSolution(std::vector<double>& vec);
			void						Rediscretization(const std::shared_ptr<Grid>&);
			void						Rediscretization();
			void						SetTimeStep(const double& step) { m_step = step; m_time = step; }
			Matrix*						GetGlobalMatrix() const;
			Grid*						GetMesh() { return m_Grid; }
			const std::vector<double>	GetRightVector() const;
			void						OutDatFormat(const Mesh::Point& min, const Mesh::Point& max, const std::string& file_name, const std::vector<double>& vec) const;
			void						OutMeshFormat(const std::string& file_name, const std::vector<double>& vec);
			void						OutMeshTimeFormat(const std::string& file_name, const std::vector<double>& vec);
			static const double	 		GetSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p);
            static const double	 		GetSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p, const int nfem);
			static const Mesh::Point	GetGradSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p);
            static const Mesh::Point	GetGradSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p, const int n);
			~FEMethod();
		private:
			void						GeneratePortrait();
			void						AssemblGlobal();
			void						MainConditions();
			void						SecondConditions();
			void						ThirdConditions();
			void						StefanConditions();
			void						ApplySources();
			const int					AssembleLocalMatrix(const int);
			Grid*						m_Grid = nullptr;
			Matrix*						m_GlobalMatrix = nullptr;
			Matrix*					m_RightMatrix = nullptr;
			Problem*					m_problem = nullptr;
			std::vector<double>			m_solution;
			std::vector<double>*		m_rhsvector;
			unsigned int				m_size;
			double						m_step{ 0.1 };
			double						m_time{ 0.1 };
			unsigned int				m_N;
			unsigned int				m_Ns;
			std::vector<unsigned int>	m_nums;

		};

		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::Discretization()
		{
			GeneratePortrait();
			AssemblGlobal();
            //ApplySources();
			//SecondConditions();
			//ThirdConditions();
            MainConditions();
			//StefanConditions();
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::GeneratePortrait()
		{
			const auto& el = m_Grid->GetElement(0);
			int order = m_Grid->GetElement(0)->GetDoFs();
			std::vector<std::set<unsigned int>> temp;
			m_Ns = m_Grid->GetNumberOfNodes();
			m_N = m_Grid->GetNumberOfElements();
			temp.resize(m_Grid->GetNumberOfNodes());
			unsigned i, j, k;
			for (k = 0; k < m_N; ++k)
			{
				const auto& elem{ m_Grid->GetElement(k) };
				for (i = 0; i < order; ++i)
					for (j = 0; j < order; ++j)
						if (elem->GetNode(j) > elem->GetNode(i))
							temp[elem->GetNode(j)].insert(elem->GetNode(i));
			}
			if(m_problem->findTerm(Terms::RUV))
				m_RightMatrix->Create(temp.size(), temp);
            
			//m_GlobalMatrix = std::shared_ptr<Matrix>(new Matrix(m_Grid->GetNumberOfNodes(), temp));
			//m_rhsvector.resize(m_Grid->GetNumberOfNodes());
            //std::cout << temp.size() << std::endl;
			m_GlobalMatrix->Create(temp.size(), temp);
			m_rhsvector->resize(temp.size());
			//m_solution.resize(m_Grid->GetNumberOfNodes());
			//for (int l = 0; l < m_Grid->GetNumberOfNodes(); ++l)
			//	m_solution[l] = 20;
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::AssemblGlobal()
		{
			int l;
			//std::vector<std::future<int>> futures;
			for (l = 0; l < m_N; ++l)
				//futures.push_back(async(&DGMethod<Problem, Grid, Matrix>::AssembleLocalMatrix, this, l));
				AssembleLocalMatrix(l);
			//for (auto &it : futures)
			//it.get();
		}
		template<class Problem, class Grid, class Matrix>
		const int FEMethod<Problem, Grid, Matrix>::AssembleLocalMatrix(const int l)
		{
			int i, j, k, nodes;
			double mij;
			const auto& elem{ m_Grid->GetElement(l) };
			const int dofs{ (int)elem->GetDoFs() };
			const int terms{ (int)m_problem->getNumberOfTerms() };
			nodes = elem->GetNumberOfNodes();
			std::vector<Mesh::Point> points(nodes);
			for (i = 0; i < nodes; ++i)
				points[i] = m_Grid->GetNode(elem->GetNode(i));
			for (k = 0; k < terms; ++k)
			{
				switch (m_problem->getTerm(k))
				{
				case Terms::IUV:
					for (i = 0; i < (int)dofs; ++i)
					{
						for (j = 0; j < (int)dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								return m_problem->get_parameter(Terms::IUV, l, elem->GetType(), p) * elem->GetShapeFunction(i, p) * elem->GetShapeFunction(j, p);
							};
							mij = elem->Integrate(M, points);
							m_GlobalMatrix->AddElement(elem->GetNode(i), elem->GetNode(j), mij);
						}
					}
					break;
				case Terms::IDUDV:
					for (i = 0; i < (int)dofs; ++i)
					{
						for (j = 0; j < (int)dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								//auto m = elem->GetGradShapeFunction(i, p) * elem->GetGradShapeFunction(j, p);
								return m_problem->get_parameter(Terms::IDUDV, l, elem->GetType(), p) * elem->GetGradShapeFunction(i, p) * elem->GetGradShapeFunction(j, p);
							};
							//mij = m_Grid->getParameter(Parameters::DIFFUSION, l, j) * elem->Integrate(M, points);
                            mij = elem->Integrate(M, points);
							m_GlobalMatrix->AddElement(elem->GetNode(i), elem->GetNode(j), mij);
						}
					}
					break;
				case Terms::IDUV:
					for (i = 0; i < (int)dofs; ++i)
					{
						for (j = 0; j < (int)dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
                                return m_problem->get_parameter(Terms::IDUV, l, elem->GetType(), p, 0) * elem->GetShapeFunction(i, p) * elem->GetGradShapeFunction(j, p);
							};
							auto _mij = elem->Integrate(M, points);
							m_GlobalMatrix->AddElement(elem->GetNode(i), elem->GetNode(j), _mij);
						}
					}
					break;
				case Terms::IUDV:
					for (i = 0; i < dofs; ++i)
					{
						for (j = 0; j < dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								return elem->GetGradShapeFunction(i, p) * elem->GetShapeFunction(j, p);
							};
							//mij = m_CoarseGrid->getParameter(Parameters::ADVECTION, l, j) * m_flux(m_CoarseGrid->getSolution(l, j)) * elem->Integrate(M, points).x;
							mij = elem->Integrate(M, points).x;
							m_GlobalMatrix->AddElement(elem->GetNode(i), elem->GetNode(j), mij);
						}
					}
					break;
				case Terms::EUV:
					for (i = 0; i < dofs; ++i)
					{
						for (j = 0; j < dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								return elem->GetShapeFunction(i, p) * elem->GetShapeFunction(j, p);
							};
							mij = elem->Integrate(M, points);
							m_rhsvector->operator[](elem->GetNode(i)) += m_Grid->getParameter(Parameters::MASS, l, j) * m_Grid->getSolution(l, j) * mij;
							//m_rhsvector->operator[](m_nums[l] + i) += m_CoarseGrid->getParameter(Parameters::MASS, l, points[j]) * elem->GetValue(j) * mij;
						}
					}
					break;
				case Terms::EDUDV:
					for (i = 0; i < dofs; ++i)
					{
						for (j = 0; j < dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								return elem->GetGradShapeFunction(i, p) * elem->GetGradShapeFunction(j, p);
							};
							mij = elem->Integrate(M, points);
							m_rhsvector->operator[](elem->GetNode(i)) += m_Grid->getParameter(Parameters::DIFFUSION, l, j) * m_Grid->getSolution(l, j) * mij;
						}
					}
					break;
				case Terms::EDUV:
					for (i = 0; i < dofs; ++i)
					{
						for (j = 0; j < dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								return elem->GetShapeFunction(i, p) * elem->GetGradShapeFunction(j, p);
							};
							mij = elem->Integrate(M, points).x;
							m_rhsvector->operator[](elem->GetNode(i)) += m_Grid->getParameter(Parameters::ADVECTION, l, j) * mij;
						}
					}
					break;
				case Terms::EUDV:
					for (i = 0; i < dofs; ++i)
					{
						for (j = 0; j < dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
								return elem->GetGradShapeFunction(i, p) * elem->GetShapeFunction(j, p);
							};
							mij = elem->Integrate(M, points).x;
							m_rhsvector->operator[](elem->GetNode(i)) += m_Grid->getParameter(Parameters::ADVECTION, l, j) * mij;// *mij;
						}
					}
					break;
				case Terms::EFV:
                    for (i = 0; i < dofs; ++i)
					{
                        /*for (j = 0; j < dofs; ++j)
                        {
                            auto M = [&](const Mesh::Point& p)
                            {
                                return m_problem->get_parameter(Terms::EFV, elem->GetType(), l, j, p) * elem->GetShapeFunction(i, p) * elem->GetShapeFunction(j, p);
                            };
                            mij = elem->Integrate(M, points);
                            m_rhsvector->operator[](elem->GetNode(i)) += mij;
                        }*/
                        auto M = [&](const Mesh::Point& p)
						{
							return m_problem->get_parameter(Terms::EFV, elem->GetType(), l, i, p) * elem->GetShapeFunction(i, p);
                        };
                        mij = elem->Integrate(M, points);
                        m_rhsvector->operator[](elem->GetNode(i)) += mij;
					}
					break;
				case Terms::RUV:
					for (i = 0; i < (int)dofs; ++i)
					{
						for (j = 0; j < (int)dofs; ++j)
						{
							auto M = [&](const Mesh::Point& p)
							{
                                return elem->GetShapeFunction(i, p) * elem->GetShapeFunction(j, p);
							};
							mij = elem->Integrate(M, points);
							m_RightMatrix->AddElement(elem->GetNode(i), elem->GetNode(j), mij);
						}
					}
                    break;
				default:
					break;
				}
			}
			return 0;
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::MainConditions()
		{
			double mu{ 1e10 };
			const auto n = m_problem->get_number_of_boundaries();
			const auto m = m_Grid->GetNumberOfBoundaries();
			for (int i = 0; i < n; ++i)
			{
				const auto& type = m_problem->get_boundary_type(i);
				for (int j = 0; j < m; ++j)
				{
					const auto& row = m_Grid->GetBoundary(j);
					if (row->GetType() == type)
					{
						const int dofs = (int)row->GetDoFs();
						const int dofs2 = 2;
						const auto& elem_num = row->GetNeighbour(0);
						const auto& elem = m_Grid->GetElement(elem_num);
						const int dofs_elem = elem->GetDoFs();
						std::vector<Mesh::Point> points(dofs_elem);
						for (int k = 0; k < dofs_elem; ++k)
							points[k] = m_Grid->GetNode(elem->GetNode(k));
                        for (int k = 0; k < dofs; ++k)
						{
							int l = 0;
							for (; l < dofs_elem; ++l)
							{
								if (elem->GetNode(l) == row->GetNode(k))
									break;
							}
							m_GlobalMatrix->NullRow(row->GetNode(k));
                            m_GlobalMatrix->operator()(row->GetNode(k), row->GetNode(k)) = 0;
							//m_GlobalMatrix->operator()(row->GetNode(k), row->GetNode(k)) *= mu;
							//m_rhsvector->operator[](row->GetNode(k)) = m_problem->get_boundary_parameter(0, type, m_Grid->GetNode(row->GetNode(k)));
							//m_rhsvector->operator[](row->GetNode(k)) = m_problem->get_boundary_parameter(0, type, elem_num, l, m_Grid->GetNode(row->GetNode(k)));
							m_rhsvector->operator[](row->GetNode(k)) = elem->GetWeight(l, points, [=](const Mesh::Point& p) { return m_problem->get_boundary_parameter(0, type, p); });
							if(m_problem->findTerm(Terms::RUV))
								{
									m_RightMatrix->NullRow(row->GetNode(k));
									//m_RightMatrix->operator()(row->GetNode(k), row->GetNode(k)) *= mu;
								}
						}
						/*for (int k = dofs2; k < dofs; ++k)
						{
							m_GlobalMatrix->NullRow(row->GetNode(k));
							m_rhsvector->operator[](row->GetNode(k)) = 0;
						}*/
					}
				}
			}
			/*for (auto bnd : m_Grid->GetBoundaryConditions())
			{
				if (get<0>(bnd.second) == 1)
					for (auto row : m_Grid->GetBoundary())
					{
						if (bnd.first == row->GetType())
						{
							for (int i = 0; i < row->GetDoF(); ++i)
							{
								m_GlobalMatrix->NullRow(row->GetNodes(i));
								m_rhsvector[row->GetNodes(i)] = get<1>(bnd.second)(m_Grid->GetNodes()[row->GetNodes(i)]);
							}
						}
					}
			}*/
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::SecondConditions()
		{
			double theta = 0;
			int nfem;
			Mesh::Point temp[3];
			std::vector<int> local;
			for (auto bnd : m_Grid->GetBoundaryConditions())
			{
				//if (get<0>(bnd.second) == 2)
				{
					for (auto row : m_Grid->GetBoundary())
					{
						if (bnd.first == row->GetType())
						{
							local.resize(0);
							int dofs = row->GetDoF();
							nfem = row->GetNumberOfElement(0);
							auto elem = m_Grid->GetElements()[nfem];
							//auto GetBasis = [&](int t, Point p){return elem->GetBasis(t, p); };
							for (int j = 0; j < dofs; ++j)
							{
								temp[j] = m_Grid->GetNodes()[row->GetNodes(j)];
								for (int i = 0; i < elem->GetDoF(); ++i)
								{
									if (row->GetNodes(j) == elem->GetNodes()[i])
									{
										local.push_back(i);
										break;
									}
								}
							}
							for (int i = 0; i < dofs; ++i)
							{
								for (int j = 0; j < dofs; ++j)
								{
									//theta = get<1>(bnd.second)(m_Grid->GetNodes()[row->GetNodes(i)]);
									theta = 0;
									auto GetMass = [&](const Mesh::Point& p) {return elem->GetBasis(local[j], p) * elem->GetBasis(local[i], p); };
									auto GetBBasis = [&](const Mesh::Point& p) {return row->GetBasis(j, p)*row->GetBasis(i, p); };
									//if (i < 2 || j < 2)
									m_rhsvector[row->GetNodes(i)] += theta * row->Integrate(GetMass, temp);

									//if (i < 3 || j < 3)
									//	m_rhsvector[row[i + 1]] += theta * row->Integrate(GetBBasis, temp);
								}
							}
						}
					}
				}
			}
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::StefanConditions()
		{
			double dest{ 0. }, lat{ 0 };
			int nfem;
			Mesh::Point temp[3];
			std::vector<int> local;
			for (auto bnd : m_Grid->GetBoundaryConditions())
			{
				//if (get<0>(bnd.second) == 4)
				{
					lat = 0;
					//lat = get<2>(bnd.second);
					for (auto row : m_Grid->GetBoundary())
					{
						if (bnd.first == row->GetType())
						{
							local.resize(0);
							int dofs = row->GetDoF();
							nfem = row->GetNumberOfElement(0);
							auto elem = m_Grid->GetElements()[nfem];
							//auto GetBasis = [&](int t, Point p){return elem->GetBasis(t, p); };
							for (int j = 0; j < dofs; ++j)
							{
								temp[j] = m_Grid->GetNodes()[row->GetNodes(j)];
								for (int i = 0; i < elem->GetDoF(); ++i)
								{
									if (row->GetNodes(j) == elem->GetNodes()[i])
									{
										local.push_back(i);
										break;
									}
								}
							}
							for (int i = 0; i < dofs; ++i)
							{
								for (int j = 0; j < dofs; ++j)
								{
									dest = 0;
									//dest = get<1>(bnd.second)(m_Grid->GetNodes()[row->GetNodes(i)]);
									auto GetBBasis = [&](const Mesh::Point& p) {return row->GetBasis(j, p)*row->GetBasis(i, p); };
									//if (i < 2 || j < 2)
									m_rhsvector[row->GetNodes(i)] += dest * lat * row->Integrate(GetBBasis, temp);

									//if (i < 3 || j < 3)
									//	m_rhsvector[row[i + 1]] += theta * row->Integrate(GetBBasis, temp);
								}
							}
						}
					}
				}
			}
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::ThirdConditions()
		{
			double param{ 0 }, beta{ 0 };
			int nfem;
			Mesh::Point temp[6];
			std::vector<int> local;
			auto fxy = [&](const Mesh::Point& p) {return (10 * p.y*m_time + m_time) / 10; };
			//auto fxy = [&](const Point& p){return 10 * p.y + 10 * m_time; };
			for (auto bnd : m_Grid->GetBoundaryConditions())
			{
				//if (get<0>(bnd.second) == 3)
				{

					for (auto row : m_Grid->GetBoundary())
					{
						if (bnd.first == row->GetType())
						{
							local.resize(0);
							int dofs = row->GetDoF();
							nfem = row->GetNumberOfElement(0);
							auto elem = m_Grid->GetElements()[nfem];
							//auto GetBasis = [&](int t, Point p){return elem->GetBasis(t, p); };
							auto order = elem->GetDoF();
							for (int j = 0; j < dofs; ++j)
							{
								temp[j] = m_Grid->GetNodes()[row->GetNodes(j)];
								for (int i = 0; i < order; ++i)
								{
									if (row->GetNodes(j) == elem->GetNodes()[i])
									{
										local.push_back(i);
										break;
									}
								}
							}
							double val{ 0 };
							for (int i = 0; i < dofs; ++i)
							{
								for (int j = 0; j < dofs; ++j)
								{
									param = 0;
									beta = 0;
									//beta = get<2>(bnd.second);
									//param = get<1>(bnd.second)(m_Grid->GetNodes()[row->GetNodes(i)]);
									//param = fxy(temp[j]);
									auto GetBBasis = [&](const Mesh::Point& p) {return elem->GetBasis(local[j], p)*elem->GetBasis(local[i], p); };
									//val = row->GetElement(GetBBasis, temp);
									val = row->Integrate(GetBBasis, temp);
									m_GlobalMatrix->operator()(row->GetNodes(i), row->GetNodes(j)) += beta * val;
									m_rhsvector[row->GetNodes(i)] += beta * param * val;
								}
							}
						}
					}
				}
			}
		}
		template<class Problem, class Grid, class Matrix>
		Matrix* FEMethod<Problem, Grid, Matrix>::GetGlobalMatrix() const
		{
			return m_GlobalMatrix;
		}
		template<class Problem, class Grid, class Matrix>
		const double FEMethod<Problem, Grid, Matrix>::GetValue(const Mesh::Point& p) const
		{
			if (!m_solution.size())
				return -1;
			double val = 0;
			int nfem = -1;
			nfem = m_Grid->FindElement(p);
			if (nfem == -1)
				return -1;
			auto elem = m_Grid->GetElements()[nfem];
			for (int i = 0; i < elem->GetDoF(); ++i)
				val += m_solution[elem->GetNodes()[i]] * elem->GetBasis(i, p);
			return val;
		}
		template<class Problem, class Grid, class Matrix>
		const double FEMethod<Problem, Grid, Matrix>::GetValue(const Mesh::Point& p, const std::vector<double>& vec) const
		{
			if (!vec.size())
				return -1;
			double val{ 0 };
			int nfem{ -1 };
			nfem = m_Grid->FindElement(p);
			if (nfem == -1)
				return -1;
			auto elem = m_Grid->GetElements()[nfem];
			for (int i = 0; i < elem->GetDoFs(); ++i)
				val += vec[elem->GetNode(i)] * elem->GetShapeFunction(i, p);
			return val;
		}
		template<class Problem, class Grid, class Matrix>
		const double FEMethod<Problem, Grid, Matrix>::GetValue(const Mesh::Point& p, const std::vector<double>& vec, const int num) const
		{
			if (!vec.size() || num < 0)
				return -1;
			double val{ 0 };
			auto elem = m_Grid->GetElements()[num];
			for (int i = 0; i < elem->GetDoF(); ++i)
				val += vec[elem->GetNodes()[i]] * elem->GetBasis(i, p);
			return val;
		}
		//template<class Problem, class Grid, class Matrix>
		//const Mesh::Point FEMethod<Problem, Grid, Matrix>::GetGradValue(const Mesh::Point& p, const std::vector<double>& vec) const
		//{
		//	Mesh::Point val{ 0, 0 };
		//	int nfem{ -1 };
		//	nfem = m_Grid->FindElement(p);
		//	if (nfem == -1)
		//		return val;
		//	auto elem = m_Grid->GetElements()[nfem];
		//	for (int i = 0; i < elem->GetDoF(); ++i)
		//	{
		//		val.x += vec[elem->GetNodes()[i]] * elem->GetGradBasis(i, p).x;
		//		val.y += vec[elem->GetNodes()[i]] * elem->GetGradBasis(i, p).y;
		//		val.z += vec[elem->GetNodes()[i]] * elem->GetGradBasis(i, p).z;
		//	}
		//	return val;
		//}
		template<class Problem, class Grid, class Matrix>
		const double FEMethod<Problem, Grid, Matrix>::GetEffective(const std::vector<double>& vec) const
		{
			double sum = 0;
			//std::vector<int> dofs;
			//Mesh::Point points[10];
			//for (int i = 0; i < m_Grid->GetElements().size(); ++i)
			//{
				//auto mb = [&](const Mesh::Point& b) {return GetGradValue(b, vec)*GetGradValue(b, vec); };
				//dofs.resize(0);
				//auto elem = m_Grid->GetElements()[i];
				//int order = elem->GetDoF();
				//double diff = std::get<0>(m_Grid->GetDiffusion().find(elem->GetType())->second);
				//for (int j = 0; j < order; ++j)
				//{
					//dofs.push_back(elem->GetNodes()[j]);
					//points[j] = m_Grid->GetNodes()[dofs[j]];
				//}
				//sum += diff * elem->Integrate(mb, points);
			//}
			//std::cout << "Effect (local): " << sum << std::endl;
			//std::cout << "Effect (local) sqrt: " << sqrt(sum) << std::endl;
			return sum;
		}
		//template<class Problem, class Grid, class Matrix>
		//const Mesh::Point FEMethod<Problem, Grid, Matrix>::GetLambdaGrad(const Mesh::Point& p, const std::vector<double>& vec) const
		//{
		//	Mesh::Point val{ 0, 0, 0 };
		//	//double val{ 0 };
		//	double diff{ 0 };
		//	Mesh::Point temp{ 0, 0, 0 };
		//	int nfem{ -1 };
		//	nfem = m_Grid->FindElement(p);
		//	if (nfem == -1)
		//		return val;
		//	auto elem = m_Grid->GetElements()[nfem];
		//	diff = std::get<0>(m_Grid->GetDiffusion().find(elem->GetType())->second);
		//	for (int i = 0; i < elem->GetDoF(); ++i)
		//	{
		//		//val += elem->GetGradBasis(i, p) * elem->GetGradBasis(i, p) * vec[elem->GetNodes()[i]] * vec[elem->GetNodes()[i]] * diff;
		//		//val += elem->GetBasis(i, p) * vec[elem->GetNodes()[i]] * diff;
		//		temp = elem->GetGradBasis(i, p);
		//		val.x += temp.x * vec[elem->GetNodes()[i]] * (diff);
		//		val.y += temp.y * vec[elem->GetNodes()[i]] * (diff);
		//		val.z += temp.z * vec[elem->GetNodes()[i]] * (diff);
		//	}
		//	return val;
		//}
		template<class Problem, class Grid, class Matrix>
		const std::vector<double> FEMethod<Problem, Grid, Matrix>::GetRightVector() const
		{
			return *m_rhsvector;
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::OutDatFormat(const Mesh::Point& mn, const Mesh::Point& mx, const std::string& file_name, const std::vector<double>& vec) const
		{
			std::ofstream of(file_name + "z.dat");
			std::streambuf *buf = std::cout.rdbuf();
			std::cout.rdbuf(of.rdbuf());
			std::cout << "TITLE = FE-METHOD\n";
			std::cout << "VARIABLES = \"dx1\", \"dx2\", \"u\"\n";
			std::cout << "ZONE i=51, j=51, F=POINT\n";
			double stepx = (mx.x - mn.x) / 51;
			double stepy = (mx.y - mn.y) / 51;
			for (int i = 0; i < 51; ++i)
				for (int j = 0; j < 51; ++j)
					std::cout << mn.x + j * stepx << "\t" << mn.y + stepy * i << "\t" << GetValue(Mesh::Point(mn.x + j * stepx, mn.y + i * stepy, mn.z), vec) << std::endl;
			std::cout.rdbuf(buf);
			of.close();
			of.open(file_name + "x.dat");
			buf = std::cout.rdbuf();
			std::cout.rdbuf(of.rdbuf());
			std::cout << "TITLE = FE-METHOD\n";
			std::cout << "VARIABLES = \"dx1\", \"dx2\", \"u\"\n";
			std::cout << "ZONE i=51, j=51, F=POINT\n";
			for (int i = 0; i < 51; ++i)
				for (int j = 0; j < 51; ++j)
					std::cout << mn.x + j * stepx << "\t" << mn.y + stepy * i << "\t" << GetValue(Mesh::Point(mn.z, mn.x + j * stepx, mn.y + i * stepy), vec) << std::endl;
			std::cout.rdbuf(buf);
			of.close();
			of.open(file_name + "y.dat");
			buf = std::cout.rdbuf();
			std::cout.rdbuf(of.rdbuf());
			std::cout << "TITLE = FE-METHOD\n";
			std::cout << "VARIABLES = \"dx1\", \"dx2\", \"u\"\n";
			std::cout << "ZONE i=51, j=51, F=POINT\n";
			for (int i = 0; i < 51; ++i)
				for (int j = 0; j < 51; ++j)
					std::cout << mn.x + j * stepx << "\t" << mn.y + stepy * i << "\t" << GetValue(Mesh::Point(mn.x + j * stepx, mn.z, mn.y + i * stepy), vec) << std::endl;
			std::cout.rdbuf(buf);
			of.close();
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::ApplySources()
		{
			int nfem = -1;
            auto total = m_problem->get_total_sources();
            for (int i = 0; i < total; ++i)
            {
                auto src = m_problem->get_point_source(i);
                auto point = src.get_point();
                nfem = m_Grid->FindElement(point);
                if (nfem != -1)
                {
                    auto val = src.get_value();
                    auto elem = m_Grid->GetElement(nfem);
                    for (int j = 0; j < 3; ++j)
                        m_rhsvector->operator[](elem->GetNode(j)) += val * elem->GetShapeFunction(j, point);
                }
                nfem = -1;
            }
            /*for (auto srd : m_Grid->GetDottedSources())
			{
				nfem = m_Grid->FindElement(srd.first);
				if (nfem != -1)
				{
					auto elem = m_Grid->GetElements()[nfem];
					for (int i = 0; i < elem->GetDoF(); ++i)
					{
						m_rhsvector[elem->GetNodes()[i]] += srd.second * elem->GetBasis(i, srd.first);
					}
				}
				nfem = -1;
            }*/
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::Rediscretization(const std::shared_ptr<Grid>& grid)
		{
			m_GlobalMatrix->NullMatrix();
			for (unsigned int i = 0; i < m_rhsvector->size(); ++i)
                (*m_rhsvector)[i] = 0;
			AssemblGlobal();
			//SecondConditions();
			//ApplySources();
			//StefanConditions();
			MainConditions();
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::Rediscretization()
		{
			m_time += m_step;
			m_GlobalMatrix->NullMatrix();
			for (unsigned int i = 0; i < m_rhsvector->size(); ++i)
                (*m_rhsvector)[i] = 0;
			AssemblGlobal();
			SecondConditions();
			ThirdConditions();
			StefanConditions();
			//ApplySources();
			MainConditions();
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::GetSolution(std::vector<double>& vec)
		{
			int size = vec.size();
			//Translation(vec);
			for (int i = 0; i < size; ++i)
				vec[i] = m_solution[i];
		}
		template<class Problem, class Grid, class Matrix>
		const double FEMethod<Problem, Grid, Matrix>::GetSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p)
		{
			double sum{ 0 };
			auto nfem{ g.FindElement(p) };
			if (nfem < 0)
				return 0.;
			auto elem{ g.GetElement(nfem) };
			auto dofs{ elem->GetDoFs() };
			for (auto i{ 0 }; i < dofs; ++i)
				sum += weights[elem->GetNode(i)] * elem->GetShapeFunction(i, p);
			return sum;
		}
        template<class Problem, class Grid, class Matrix>
        const double FEMethod<Problem, Grid, Matrix>::GetSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p, const int nfem)
        {
            double sum{ 0 };
            //if (nfem < 0)
            //    return 0.;
            auto elem{ g.GetElement(nfem) };
            auto dofs{ elem->GetDoFs() };
            //std::cout << nfem << std::endl;
            for (auto i{ 0 }; i < dofs; ++i)
                sum += weights[elem->GetNode(i)] * elem->GetShapeFunction(i, p);
            return sum;
        }
		template<class Problem, class Grid, class Matrix>
		const Mesh::Point FEMethod<Problem, Grid, Matrix>::GetGradSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p)
		{
			Mesh::Point sum{ 0, 0, 0 };
			auto nfem{ g.FindElement(p) };
			auto elem{ g.GetElement(nfem) };
			auto dofs{ elem->GetDoFs() };
			for (auto i{ 0 }; i < dofs; ++i)
				sum += weights[elem->GetNode(i)] * elem->GetGradShapeFunction(i, p);
			return sum;
		}
        template<class Problem, class Grid, class Matrix>
        const Mesh::Point FEMethod<Problem, Grid, Matrix>::GetGradSolution(const Grid& g, const std::vector<double> &weights, const Mesh::Point& p, const int nfem)
        {
            Mesh::Point sum{ 0, 0, 0 };
            auto elem{ g.GetElement(nfem) };
            auto dofs{ elem->GetDoFs() };
            for (auto i{ 0 }; i < dofs; ++i)
                sum += weights[elem->GetNode(i)] * elem->GetGradShapeFunction(i, p);
            return sum;
        }
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::LoadSolution(const std::vector<double>& vec)
		{
			m_solution.resize(vec.size());
			for (unsigned int i = 0; i < vec.size(); ++i)
				m_solution[i] = vec[i];
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::OutMeshFormat(const std::string& file_name, const std::vector<double>& vec)
		{
			const int size{ (int)m_Grid->GetNodes().size() };
			const int number{ (int)m_Grid->GetElements().size() };
			//const int size{ number * 4 };
			std::ofstream ofs(file_name + ".dat", std::ios::out);
			std::string title("TITLE = \"Mesh data\"\n Variables = \"X\", \"Y\", \"Z\", \"U\"\n Zone N = " + std::to_string(size) + ", E = " + std::to_string(number) + ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON\n");
			ofs << title;
			Mesh::Point p;
			for (int i = 0; i < size; ++i)
			{
				p = m_Grid->GetNodes()[i];
				ofs << p.x << "\t" << p.y << "\t" << p.z << "\t" << GetValue(p, vec, 1) << std::endl;
			}
			for (int i = 0; i < number; ++i)
			{
				auto elem = m_Grid->GetElements()[i];
				for (int k = 0; k < 4; ++k)
				{
					ofs << elem->GetNodes()[k] + 1 << "\t";
				}
				ofs << std::endl;
			}
			ofs.close();
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::OutMeshTimeFormat(const std::string& file_name, const std::vector<double>& vec)
		{
			const int size{ (int)m_Grid->GetNodes().size() };
			const int number{ (int)m_Grid->GetElements().size() };
			//const int size{ number * 4 };
			std::ofstream ofs(file_name + ".dat", std::ios::out | std::ios::app);
			std::string title("TITLE = \"Mesh data\"\n Variables = \"X\", \"Y\", \"Z\", \"U\"\n Zone N = " + std::to_string(size) + ", E = " + std::to_string(number) + ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON\n");
			ofs << title;
			Mesh::Point p;
			for (int i = 0; i < size; ++i)
			{
				p = m_Grid->GetNodes()[i];
				ofs << p.x << "\t" << p.y << "\t" << p.z << "\t" << GetValue(p, vec, 1) << std::endl;
			}
			for (int i = 0; i < number; ++i)
			{
				auto elem = m_Grid->GetElements()[i];
				for (int k = 0; k < 4; ++k)
				{
					ofs << elem->GetNodes()[k] + 1 << "\t";
				}
				ofs << std::endl;
			}
			ofs.close();
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::ProjectSolution(std::vector<double>& sol, std::function<const double(const Mesh::Point&, const std::vector<double>&, const int)> GetVal, std::vector<double>& vec)
		{
			for (int i = 0; i < m_Grid->GetElements().size(); ++i)
			{
				auto elem = m_Grid->GetElements()[i];
				int order = elem->GetDoF();
				for (int j = 0; j < order; ++j)
					sol[elem->GetNodes(j)] = GetVal(m_Grid->GetNodes()[elem->GetNodes(j)], vec, i);
			}
		}
		template<class Problem, class Grid, class Matrix>
		void FEMethod<Problem, Grid, Matrix>::ProjectSolution(std::vector<double>& sol, std::function<const double(const Mesh::Point&, const std::vector<double>&)> GetVal, std::vector<double>& vec, const int)
		{
			for (int i = 0; i < m_Grid->GetElements().size(); ++i)
			{
				auto elem = m_Grid->GetElements()[i];
				int order = elem->GetDoF();
				for (int j = 0; j < order; ++j)
					sol[elem->GetNodes(j)] = GetVal(m_Grid->GetNodes()[elem->GetNodes(j)], vec);
			}
		}
		template<class Problem, class Grid, class Matrix>
		const std::vector<double> FEMethod<Problem, Grid, Matrix>::SetSolution(const int sol, const int liq, const double s, const double l, const double m)
		{
			int i;
			m_solution.resize(m_Grid->GetNodes().size());
			for (i = 0; i < m_Grid->GetElements().size(); ++i)
			{
				auto elem = m_Grid->GetElements()[i];
				int order = elem->GetDoF();
				if (m_Grid->GetElements()[i]->GetType() == liq)
					for (int j = 0; j < order; ++j)
						m_solution[elem->GetNodes()[j]] = l;
				else
					for (int j = 0; j < order; ++j)
						m_solution[elem->GetNodes()[j]] = s;
			}

			for (auto bnd : m_Grid->GetBoundaryConditions())
			{
				//if (get<0>(bnd.second) == 4)
				{
					for (auto row : m_Grid->GetBoundary())
					{
						if (bnd.first == row->GetType())
						{
							int dofs = row->GetDoF();
							for (int i = 0; i < dofs; ++i)
							{
								m_solution[row->GetNodes(i)] = m;
							}
						}
					}
				}
			}
			return m_solution;
		}
		template<class Problem, class Grid, class Matrix>
		FEMethod<Problem, Grid, Matrix>::~FEMethod()
		{
			delete m_Grid;
		}
	}
}

#endif // !CORENC_METHODS_FEMethod_h

