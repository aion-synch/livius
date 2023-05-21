#ifndef CORENC_SOLVERS_DG_SOLVER_SHALLOW_WATER_H_
#define CORENC_SOLVERS_DG_SOLVER_SHALLOW_WATER_H_

#include <vector>
#include <functional>
#include <istream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../CoreNCFEM/Point.h"
namespace corenc
{
	namespace solvers
	{
		struct vector_solution
		{
			std::vector<double> S[3];
			vector_solution() {}
			vector_solution(const int _size)
			{
				S[0].resize(_size);
				S[1].resize(_size);
				S[2].resize(_size);
			}
		};
		class dg_solver_shallow_water
		{
		public:
			dg_solver_shallow_water();
			~dg_solver_shallow_water();
			const int				solve() const;
			const int				solve(
				const double t0,
				const double t1,
				const size_t nx,
				const size_t ny,
				const double x0,
				const double x1,
				const double y0,
				const double y1,
				const double g,
				const double H,
				const std::function<const std::vector<double>(const std::vector<double>&)>&, 
				const std::function<const std::vector<double>(const std::vector<double>&)>&,
				const std::function<const std::vector<double>(const std::vector<double>&)>&) const;
		};
		
		template<class Mesh>
		class dg_shallow_water
		{
		public:
			dg_shallow_water();
			~dg_shallow_water();
			const int				solve(
				const double t0,
				const double t1,
				const Mesh& mesh,
				vector_solution& sol,
				const std::function<const std::vector<double>(const std::vector<double>&)>&,
				const std::function<const std::vector<double>(const std::vector<double>&)>&,
				const std::function<const std::vector<double>(const std::vector<double>&)>&) const;
			const int				solve(
				const double t0,
				const double t1,
				const Mesh& mesh,
				vector_solution& sol,
				std::vector<double>& bath,
				std::vector<double>& ze,
				std::vector<double>& dzx,
				std::vector<double>& dzy,
				std::vector<double>& dbx,
				std::vector<double>& dby,
				const std::function<const std::vector<double>(const std::vector<double>&, const int)>&,
				const std::function<const std::vector<double>(const std::vector<double>&, const int)>&,
				const std::function<const std::vector<double>(const std::vector<double>&, const int)>&,
				const bool WRITE_FILE) const;
		};
		
		template<class Mesh>
		dg_shallow_water<Mesh>::dg_shallow_water()
		{

		}
		template<class Mesh>
		dg_shallow_water<Mesh>::~dg_shallow_water()
		{

		}
		
		template<class Mesh>
		const int dg_shallow_water<Mesh>::solve(
			const double t0,
			const double t1,
			const Mesh& mesh,
			vector_solution& sol,
			const std::function < const std::vector<double>(const std::vector<double>&)>&R,
			const std::function < const std::vector<double>(const std::vector<double>&)>&G,
			const std::function < const std::vector<double>(const std::vector<double>&)>&F) const
		{
			std::vector<double> Ut[3];
			const int max_iter = 30000;
			const double dx = mesh.GetNode(mesh.GetNumberOfNodes() - 1).x - mesh.GetNode(0).x;
			const double dy = mesh.GetNode(mesh.GetNumberOfNodes() - 1).y - mesh.GetNode(0).y;
			//const double dx = (x1 - x0) / nx;
			//const double dy = (y1 - y0) / ny;
			const int size = mesh.GetNumberOfElements();
			const int bsize = mesh.GetNumberOfBoundaries();
			
			std::vector<vector_solution> U(2);
			std::vector<vector_solution> W(2);
			U[0].S[0].resize(size);
			U[0].S[1].resize(size);
			U[0].S[2].resize(size);
			U[1].S[0].resize(size);
			U[1].S[1].resize(size);
			U[1].S[2].resize(size);

			W[0].S[0].resize(size);
			W[0].S[1].resize(size);
			W[0].S[2].resize(size);
			W[1].S[0].resize(size);
			W[1].S[1].resize(size);
			W[1].S[2].resize(size);
			for (size_t i = 0; i < size; ++i)
			{
				W[0].S[0][i] = sol.S[0][i];
				W[0].S[1][i] = sol.S[1][i];
				W[0].S[2][i] = sol.S[2][i];

				U[0].S[0][i] = sol.S[0][i];
				U[0].S[1][i] = sol.S[1][i] / sol.S[0][i];
				U[0].S[2][i] = sol.S[2][i] / sol.S[0][i];
			}

			double t_step = 0.1;
			const double cfl = 0.5;
			// W = [h hu hv]
			double lambda_x = 0;
			double lambda_y = 0;
			double lambdax = 0;
			double lambday = 0;
			double lambda = 0;
			double t_curr = 0;
			double g = 1;
			size_t iter_max = 10000;
			for (size_t t = 0; t < iter_max && t_curr < t1; ++t, t_curr += t_step)
			{
				lambda_x = 0;
				lambda_y = 0;
				for (size_t i = 0; i < size; ++i)
				{
					const auto& elem = mesh.GetElement(i);
					const auto& res = F(std::vector<double>{W[t].S[0][i], W[t].S[1][i], W[t].S[2][i]});
					W[t + 1].S[0][i] = W[t].S[0][i] + res[0];
					W[t + 1].S[1][i] = W[t].S[1][i] + res[1];
					W[t + 1].S[2][i] = W[t].S[2][i] + res[2];

					lambda_x = std::max(fabs(U[t].S[1][i]), lambda_x);
					lambda_y = std::max(fabs(U[t].S[2][i]), lambda_y);
				}
				t_step = cfl / 2 * std::min(dx / lambda_x, dy / lambda_y);
				if (t_curr + t_step > t1)
					t_step = t1 - t_curr;
				for (size_t i = 0; i < bsize; ++i)
				{
					const auto& bound = mesh.GetBoundary(i);
					const int nk = bound->GetNeighbour(0);
					const int ne = bound->GetNeighbour(1);
					const auto& normal = bound->GetNormal();
					if (ne > -1)
					{
						const auto& normal = bound->GetNormal();
						std::vector<double> wk(3);
						wk[0] = U[t].S[0][nk];
						wk[1] = U[t].S[1][nk] * U[t].S[0][nk];
						wk[2] = U[t].S[2][nk] * U[t].S[0][nk];

						std::vector<double> we(3);
						we[0] = U[t].S[0][ne];
						we[1] = U[t].S[1][ne] * U[t].S[0][ne];
						we[2] = U[t].S[2][ne] * U[t].S[0][ne];

						//lambda_x = std::max(fabs(U[t].S[1][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[1][ne]) + sqrt(g * U[t].S[0][ne]));
						//lambda_y = std::max(fabs(U[t].S[2][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[2][ne]) + sqrt(g * U[t].S[0][ne]));

						lambda_x = std::max(fabs(U[t].S[1][nk]), fabs(U[t].S[1][ne]));
						lambda_y = std::max(fabs(U[t].S[2][nk]), fabs(U[t].S[2][ne]));


						lambdax = std::max(lambdax, lambda_x);
						lambday = std::max(lambday, lambda_y);
						double ll = std::max(lambda_x, lambda_y);
						//cout << "max:\t" << ll << endl;
						std::vector<double> uk(3);
						uk[0] = U[t].S[0][nk];
						uk[1] = U[t].S[1][nk];
						uk[2] = U[t].S[2][nk];

						std::vector<double> ue(3);
						ue[0] = U[t].S[0][ne];
						ue[1] = U[t].S[1][ne];
						ue[2] = U[t].S[2][ne];

						const auto rk = R(uk);
						const auto re = R(ue);
						const auto gk = G(uk);
						const auto ge = G(ue);

						std::vector<double> uu(3);
						uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
						uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
						uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						W[t + 1].S[0][nk] -= uu[0];
						W[t + 1].S[1][nk] -= uu[1];
						W[t + 1].S[2][nk] -= uu[2];

						uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[0] + re[0]) / 2 - normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
						uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[1] + re[1]) / 2 - normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
						uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[2] + re[2]) / 2 - normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						W[t + 1].S[0][ne] -= uu[0];
						W[t + 1].S[1][ne] -= uu[1];
						W[t + 1].S[2][ne] -= uu[2];
					}
				}
				for (size_t i = 0; i < bsize; ++i)
				{
					const auto& bound = mesh.GetBoundary(i);
					const int nk = bound->GetNeighbour(0);
					const int ne = bound->GetNeighbour(1);
					if (ne == -1)
					{
						auto normal = bound->GetNormal();

						std::vector<double> u(3);
						u[0] = U[t].S[0][nk];
						u[1] = -U[t].S[1][nk];
						u[2] = -U[t].S[2][nk];

						lambdax = std::max(lambdax, lambda_x);
						lambday = std::max(lambday, lambda_y);

						std::vector<double> uk(3);
						uk[0] = U[t].S[0][nk];
						uk[1] = U[t].S[1][nk];
						uk[2] = U[t].S[2][nk];

						std::vector<double> ue(3);
						ue[0] = u[0];
						ue[1] = u[1];
						ue[2] = u[2];

						const auto rk = R(uk);
						const auto re = R(ue);
						const auto gk = G(uk);
						const auto ge = G(ue);
						lambda_x = std::max(fabs(uk[1]), fabs(ue[1]));
						lambda_y = std::max(fabs(uk[2]), fabs(ue[2]));
						std::vector<double> uu(3);
						if (normal.x > 0 || normal.y > 0)
						{
							uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
							uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
							uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						}
						else
						{
							uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
							uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
							uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						}

						W[t + 1].S[0][nk] -= uu[0];
						W[t + 1].S[1][nk] -= uu[1];
						W[t + 1].S[2][nk] -= uu[2];
					}
				}
				for (size_t i = 0; i < size; ++i)
				{

					U[t + 1].S[0][i] = W[t + 1].S[0][i];
					U[t + 1].S[1][i] = W[t + 1].S[1][i] / W[t + 1].S[0][i];
					U[t + 1].S[2][i] = W[t + 1].S[2][i] / W[t + 1].S[0][i];
				}
				W.push_back(vector_solution(size));
				U.push_back(vector_solution(size));

			}
			const auto ut = W.size() - 2;
			for (size_t i = 0; i < size; ++i)
			{
				sol.S[0][i] = W[ut].S[0][i];
				sol.S[1][i] = W[ut].S[1][i];
				sol.S[2][i] = W[ut].S[2][i];
			}
			return 0;
		}

		template<class Mesh>
		const int dg_shallow_water<Mesh>::solve(
			const double t0,
			const double t1,
			const Mesh& mesh,
			vector_solution& sol,
			std::vector<double>& bath,
			std::vector<double>& ze,
			std::vector<double>& dzx,
			std::vector<double>& dzy,
			std::vector<double>& dbx,
			std::vector<double>& dby,
			const std::function < const std::vector<double>(const std::vector<double>&, const int)>&R,
			const std::function < const std::vector<double>(const std::vector<double>&, const int)>&G,
			const std::function < const std::vector<double>(const std::vector<double>&, const int)>&F,
			const bool WRITE_FILE) const
		{
			std::vector<double> Ut[3];
			const int max_iter = 30000;
			double dx = 100, dy = 100;
			//const double dx = mesh.GetNode(mesh.GetNumberOfNodes() - 1).x - mesh.GetNode(0).x;
			//const double dy = mesh.GetNode(mesh.GetNumberOfNodes() - 1).y - mesh.GetNode(0).y;
			//const double dx = (x1 - x0) / nx;
			//const double dy = (y1 - y0) / ny;
			const int size = mesh.GetNumberOfElements();
			const int bsize = mesh.GetNumberOfBoundaries();

			std::vector<vector_solution> U(2);
			std::vector<vector_solution> W(2);
			U[0].S[0].resize(size);
			U[0].S[1].resize(size);
			U[0].S[2].resize(size);
			U[1].S[0].resize(size);
			U[1].S[1].resize(size);
			U[1].S[2].resize(size);

			W[0].S[0].resize(size);
			W[0].S[1].resize(size);
			W[0].S[2].resize(size);
			W[1].S[0].resize(size);
			W[1].S[1].resize(size);
			W[1].S[2].resize(size);
			for (size_t i = 0; i < size; ++i)
			{
				W[0].S[0][i] = sol.S[0][i];
				W[0].S[1][i] = sol.S[1][i];
				W[0].S[2][i] = sol.S[2][i];

				U[0].S[0][i] = sol.S[0][i];
				U[0].S[1][i] = sol.S[1][i] / sol.S[0][i];
				U[0].S[2][i] = sol.S[2][i] / sol.S[0][i];
			}
			auto center = [=](const size_t i)
			{
				const auto& elem = mesh.GetElement(i);
				std::vector<corenc::Mesh::Point> pts(4);
				pts[0] = mesh.GetNode(elem->GetNode(0));
				pts[1] = mesh.GetNode(elem->GetNode(1));
				pts[2] = mesh.GetNode(elem->GetNode(2));
				pts[3] = mesh.GetNode(elem->GetNode(3));
				return corenc::Mesh::Point(pts[0].x + (pts[3].x - pts[0].x) / 2, pts[0].y + (pts[3].y - pts[0].y) / 2);
			};
			double t_step = 0.1;
			const double cfl = 0.1;
			// W = [h hu hv]
			double lambda_x = 0;
			double lambda_y = 0;
			double lambdax = 0;
			double lambday = 0;
			double lambda = 0;
			double t_curr = 0;
			double g = 1;
			size_t iter_max = 100;
			for (size_t t = 0; t < iter_max && t_curr < t1; ++t, t_curr += t_step)
			{
				lambda_x = 0;
				lambda_y = 0;
				for (size_t i = 0; i < size; ++i)
				{
					const auto& elem = mesh.GetElement(i);
					const auto& res = F(std::vector<double>{W[t].S[0][i], W[t].S[1][i], W[t].S[2][i]}, i);
					W[t + 1].S[0][i] = W[t].S[0][i] + res[0];
					W[t + 1].S[1][i] = W[t].S[1][i] + res[1];
					W[t + 1].S[2][i] = W[t].S[2][i] + res[2];

					lambda_x = std::max(fabs(U[t].S[1][i]) + sqrt(g*U[t].S[0][i]), lambda_x);
					lambda_y = std::max(fabs(U[t].S[2][i]) + sqrt(g*U[t].S[0][i]), lambda_y);
					dx = std::min(mesh.GetNode(elem->GetNode(3)).x - mesh.GetNode(elem->GetNode(0)).x, dx);
					dy = std::min(mesh.GetNode(elem->GetNode(3)).y - mesh.GetNode(elem->GetNode(0)).y, dy);
					//lambda_x = std::min(U[t].S[1][i])
					//lambda_x = std::max(fabs(U[t].S[1][i]), lambda_x);
					//lambda_y = std::max(fabs(U[t].S[2][i]), lambda_y);
				}
				t_step = cfl / 2 * std::min(dx / lambda_x, dy / lambda_y);
				//std::cout << t_step << std::endl;
				if (t_curr + t_step > t1)
					t_step = t1 - t_curr;
				for (size_t i = 0; i < bsize; ++i)
				{
					const auto& bound = mesh.GetBoundary(i);
					const int nk = bound->GetNeighbour(0);
					const int ne = bound->GetNeighbour(1);
					const auto& normal = bound->GetNormal();
					if (ne > -1)
					{
						const auto& normal = bound->GetNormal();
						std::vector<double> wk(3);
						wk[0] = U[t].S[0][nk];
						wk[1] = U[t].S[1][nk] * U[t].S[0][nk];
						wk[2] = U[t].S[2][nk] * U[t].S[0][nk];

						std::vector<double> we(3);
						we[0] = U[t].S[0][ne];
						we[1] = U[t].S[1][ne] * U[t].S[0][ne];
						we[2] = U[t].S[2][ne] * U[t].S[0][ne];

						lambda_x = std::max(fabs(U[t].S[1][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[1][ne]) + sqrt(g * U[t].S[0][ne]));
						lambda_y = std::max(fabs(U[t].S[2][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[2][ne]) + sqrt(g * U[t].S[0][ne]));

						//lambda_x = std::max(fabs(U[t].S[1][nk]), fabs(U[t].S[1][ne]));
						//lambda_y = std::max(fabs(U[t].S[2][nk]), fabs(U[t].S[2][ne]));


						lambdax = std::max(lambdax, lambda_x);
						lambday = std::max(lambday, lambda_y);
						double ll = std::max(lambda_x, lambda_y);
						//cout << "max:\t" << ll << endl;
						std::vector<double> uk(3);
						uk[0] = U[t].S[0][nk];
						uk[1] = U[t].S[1][nk];
						uk[2] = U[t].S[2][nk];

						std::vector<double> ue(3);
						ue[0] = U[t].S[0][ne];
						ue[1] = U[t].S[1][ne];
						ue[2] = U[t].S[2][ne];

						const auto rk = R(uk, nk);
						const auto re = R(ue, ne);
						const auto gk = G(uk, nk);
						const auto ge = G(ue, ne);

						std::vector<double> uu(3);
						uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
						uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
						uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						W[t + 1].S[0][nk] -= uu[0];
						W[t + 1].S[1][nk] -= uu[1];
						W[t + 1].S[2][nk] -= uu[2];

						uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[0] + re[0]) / 2 - normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
						uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[1] + re[1]) / 2 - normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
						uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[2] + re[2]) / 2 - normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						W[t + 1].S[0][ne] -= uu[0];
						W[t + 1].S[1][ne] -= uu[1];
						W[t + 1].S[2][ne] -= uu[2];


					}
				}
				for (size_t i = 0; i < bsize; ++i)
				{
					const auto& bound = mesh.GetBoundary(i);
					const int nk = bound->GetNeighbour(0);
					const int ne = bound->GetNeighbour(1);
					if (ne == -1)
					{
						auto normal = bound->GetNormal();

						std::vector<double> u(3);
						u[0] = U[t].S[0][nk];
						u[1] = -U[t].S[1][nk];
						u[2] = -U[t].S[2][nk];

						//u[0] = U[t].S[0][nk];
						//u[1] = U[t].S[1][nk];
						//u[2] = U[t].S[2][nk];

						lambdax = std::max(lambdax, lambda_x);
						lambday = std::max(lambday, lambda_y);

						std::vector<double> uk(3);
						uk[0] = U[t].S[0][nk];
						uk[1] = U[t].S[1][nk];
						uk[2] = U[t].S[2][nk];

						std::vector<double> ue(3);
						ue[0] = u[0];
						ue[1] = u[1];
						ue[2] = u[2];

						const auto rk = R(uk, nk);
						const auto re = R(ue, nk);
						const auto gk = G(uk, nk);
						const auto ge = G(ue, nk);
						lambda_x = std::max(fabs(uk[1]) + sqrt(g*uk[0]), fabs(ue[1]) + sqrt(g*ue[0]));
						lambda_y = std::max(fabs(uk[2]) + sqrt(g*uk[0]), fabs(ue[2]) + sqrt(g*ue[0]));
						std::vector<double> uu(3);
						if (normal.x > 0 || normal.y > 0)
						{
							uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
							uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
							uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						}
						else
						{
							uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
							uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
							uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
						}

						W[t + 1].S[0][nk] -= uu[0];
						W[t + 1].S[1][nk] -= uu[1];
						W[t + 1].S[2][nk] -= uu[2];
					}
				}
				for (size_t i = 0; i < size; ++i)
				{

					U[t + 1].S[0][i] = W[t + 1].S[0][i];
					U[t + 1].S[1][i] = W[t + 1].S[1][i] / W[t + 1].S[0][i];
					U[t + 1].S[2][i] = W[t + 1].S[2][i] / W[t + 1].S[0][i];
				}
				W.push_back(vector_solution(size));
				U.push_back(vector_solution(size));

				for (size_t k = 0; k < bsize; ++k)
				{
					const auto& bound = mesh.GetBoundary(k);
					const auto nk = bound->GetNeighbour(0);
					const auto ne = bound->GetNeighbour(1);
					ze[nk] = W[t + 1].S[0][nk] - bath[nk];
					if (ne > -1)
					{
						ze[ne] = W[t + 1].S[0][ne] - bath[ne];
						const auto ce = center(ne);
						const auto ck = center(nk);
						const double cx = ce.x - ck.x;
						const double cy = ce.y - ck.y;
						if (fabs(cy) < 1e-13)
						{
							dzx[nk] = (ze[ne] - ze[nk]) / cx;
							dzx[ne] = dzx[nk];
							dbx[nk] = (bath[ne] - bath[nk]) / cx;
							dbx[ne] = dbx[nk];
						}
						else
						{
							dzy[nk] = (ze[ne] - ze[nk]) / cy;
							dzy[ne] = dzy[nk];
							dby[nk] = (bath[ne] - bath[nk]) / cy;
							dby[ne] = dby[nk];
						}
					}
				}
			}
			/*const auto ut = W.size() - 2;
			for (size_t i = 0; i < size; ++i)
			{
				sol.S[0][i] = W[ut].S[0][i];
				sol.S[1][i] = W[ut].S[1][i];
				sol.S[2][i] = W[ut].S[2][i];
			}
			std::ofstream ofs;
			ofs.open("meshU.txt");
			const size_t t_r = U.size() - 1;
			ofs << t_r << std::endl;
			for (size_t i = 0; i < t_r; ++i)
				for (size_t j = 0; j < size; ++j)
					ofs << U[i].S[0][j] - bath[j] << std::endl;
			ofs.close();*/
			return 0;
		}
	}
}


#endif // !CORENC_SOLVERS_DG_SOLVER_SHALLOW_WATER_H_
