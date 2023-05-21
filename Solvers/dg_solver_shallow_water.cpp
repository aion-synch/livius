#include "dg_solver_shallow_water.h"
#include <vector>
#include "../CoreNCFEM/Grids/RegularMesh.h"
#include "../CoreNCFEM/Parameter.h"
#include <algorithm>
#include <functional>
using namespace corenc;
using namespace std;
using namespace solvers;
using namespace Mesh;
dg_solver_shallow_water::dg_solver_shallow_water()
{

}

dg_solver_shallow_water::~dg_solver_shallow_water()
{

}
const int dg_solver_shallow_water::solve() const
{
	// U = [H Hu Hv]
	// R = [Hu Hu^2 Huv]
	// G = [Hv Huv Hv^2]
	// S = [0 -gHdkdx+1/p(t_sx-t_bx+Fx -gHdkdy+1/p(t_sy-t_by+Fy))]
	// test
	// U = [h hu hv]
	// R = [hu hu^2+gh^2/2 huv]
	// G = [hv huv hv^2+gh^2/2]
	// S = [0 0 0]
	vector<double> Ut[3];
	const double t0 = 0;
	const double t1 = 3;
	const int max_iter = 30000;
	const size_t nx = 4;
	const size_t ny = 4;
	const double x0 = -1;
	const double x1 = 1;
	const double y0 = -1;
	const double y1 = 1;
	const double dx = (x1 - x0) / nx;
	const double dy = (y1 - y0) / ny;
	CRegularMesh mesh{ x0, y0, x1, y1, nx, ny };
	// dunno
	auto center = [=](const int i)
	{
		const auto& elem = mesh.GetElement(i);
		vector<Point> pts(4);
		pts[0] = mesh.GetNode(elem->GetNode(0));
		pts[1] = mesh.GetNode(elem->GetNode(1));
		pts[2] = mesh.GetNode(elem->GetNode(2));
		pts[3] = mesh.GetNode(elem->GetNode(3));
		return Point(pts[0].x + (pts[3].x - pts[0].x) / 2, pts[0].y + (pts[3].y - pts[0].y) / 2);
	};
	const int size = mesh.GetNumberOfElements();
	const int bsize = mesh.GetNumberOfBoundaries();
	struct solution
	{
		vector<double> S[3];
		solution(){}
		solution(const int _size) 
		{ 
			S[0].resize(_size); 
			S[1].resize(_size);
			S[2].resize(_size);
		}
	};
	vector<solution> U(2);
	vector<solution> W(2);
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
		const auto cent(center(i));
		if (cent.x >= -1. / 2 && cent.x <= 1. / 2 && cent.y >= -1. / 2 && cent.y <= 1. / 2)
		{
			U[0].S[0][i] = 4;
			U[0].S[1][i] = 0;
			U[0].S[2][i] = 0;

			W[0].S[0][i] = 4;
			W[0].S[1][i] = 0;
			W[0].S[2][i] = 0;
		}
		else
		{
			U[0].S[0][i] = 1;
			U[0].S[1][i] = 0;
			U[0].S[2][i] = 0;

			W[0].S[0][i] = 1;
			W[0].S[1][i] = 0;
			W[0].S[2][i] = 0;
		}
	}
	//parameter<double> g(1);
	double g = 1;
	double H = 1;
	//parameter<double> H(0);
	parameter<double> p(1);
	parameter<double> tau_sx(0);
	parameter<double> tau_sy(0);
	parameter<double> tau_bx(0);
	parameter<double> tau_by(0);
	parameter<double> Fx(0);
	parameter<double> Fy(0);
	double t_step = 0.1;
	const double cfl = 0.8;
	// W = [h hu hv]
	auto R = [=](const vector<double>& _u)
	{
		// U = [h u v]
		// R = [hu hu^2+gh^2/2 huv]
		vector<double> r(3);
		const double h = _u[0];
		const double u = _u[1];
		const double v = _u[2];
		r[0] = h * u;
		r[1] = h * u * u + g * h * h / 2;
		r[2] = h * u * v;
		return r;
	};
	auto G = [=](const vector<double>& _u)
	{
		// U = [h u v]
		// G = [hv huv hv^2+gh^2/2]
		vector<double> r(3);
		const double h = _u[0];
		const double u = _u[1];
		const double v = _u[2];
		r[0] = h * v;
		r[1] = h * u * v;
		r[2] = h * v * v + g * h * h / 2;
		return r;
	};
	double lambda_x = 0;
	double lambda_y = 0;
	double lambdax = 0;
	double lambday = 0;
	double t_curr = 0;
	size_t iter_max = 10000;
	for (size_t t = 0; t < iter_max && t_curr < t1; ++t, t_curr += t_step)
	{
		cout << "dt:\t" << t_step << endl;
		cout << "t:\t" << t_curr << endl;
		lambda_x = 0;
		lambda_y = 0;
		for (size_t i = 0; i < size; ++i)
		{
			W[t + 1].S[0][i] = W[t].S[0][i];
			W[t + 1].S[1][i] = W[t].S[1][i];// *U[t].S[1][i];
			W[t + 1].S[2][i] = W[t].S[2][i];// *U[t].S[2][i];

			//U[t + 1].S[0][i] = U[t].S[0][i];
			//U[t + 1].S[1][i] = U[t].S[1][i];
			//U[t + 1].S[2][i] = U[t].S[2][i];

			lambda_x = std::max(fabs(U[t].S[1][i]) + sqrt(g * U[t].S[0][i]), lambda_x);
			lambda_y = std::max(fabs(U[t].S[2][i]) + sqrt(g * U[t].S[0][i]), lambda_y);
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
				vector<double> wk(3);
				wk[0] = U[t].S[0][nk];
				wk[1] = U[t].S[1][nk] * U[t].S[0][nk];
				wk[2] = U[t].S[2][nk] * U[t].S[0][nk];

				vector<double> we(3);
				we[0] = U[t].S[0][ne];
				we[1] = U[t].S[1][ne] * U[t].S[0][ne];
				we[2] = U[t].S[2][ne] * U[t].S[0][ne];

				//lambda_x = U[t].S[1][nk] / U[t].S[1][nk] + sqrt(g * U[t].S[0][nk]);
				//lambda_y = U[t].S[2][nk] / U[t].S[0][nk] + sqrt(g * U[t].S[0][ne]);

				lambda_x = std::max(fabs(U[t].S[1][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[1][ne]) + sqrt(g * U[t].S[0][ne]));
				lambda_y = std::max(fabs(U[t].S[2][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[2][ne]) + sqrt(g * U[t].S[0][ne]));

				lambdax = std::max(lambdax, lambda_x);
				lambday = std::max(lambday, lambda_y);
				double ll = max(lambda_x, lambda_y);
				//cout << "max:\t" << ll << endl;
				vector<double> uk(3);
				uk[0] = U[t].S[0][nk];
				uk[1] = U[t].S[1][nk];
				uk[2] = U[t].S[2][nk];

				vector<double> ue(3);
				ue[0] = U[t].S[0][ne];
				ue[1] = U[t].S[1][ne];
				ue[2] = U[t].S[2][ne];

				const auto rk = R(uk);
				const auto re = R(ue);
				const auto gk = G(uk);
				const auto ge = G(ue);

				vector<double> uu(3);
				uu[0] = t_step / dx * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
				uu[1] = t_step / dx * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
				uu[2] = t_step / dx * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
				if (nk == 7)
				{
					//cout << "nk" << endl;
					//cout << uu[1] << endl;
					//cout << uu[1] << endl;
					//cout << uu[2] << endl;
				}
				W[t + 1].S[0][nk] -= uu[0];
				W[t + 1].S[1][nk] -= uu[1];
				W[t + 1].S[2][nk] -= uu[2];

				//U[t + 1].S[0][nk] -= uu[0];
				//U[t + 1].S[1][nk] -= uu[1] / uu[0];
				//U[t + 1].S[2][nk] -= uu[2] / uu[0];

				uu[0] = t_step / dx * (-normal.x * (rk[0] + re[0]) / 2 - normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
				uu[1] = t_step / dx * (-normal.x * (rk[1] + re[1]) / 2 - normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
				uu[2] = t_step / dx * (-normal.x * (rk[2] + re[2]) / 2 - normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
				W[t + 1].S[0][ne] -= uu[0];
				W[t + 1].S[1][ne] -= uu[1];
				W[t + 1].S[2][ne] -= uu[2];
				if (ne == 7)
				{
					//cout << "ne" << endl;
					//cout << uu[1] << endl;
					//cout << uu[1] << endl;
					//cout << uu[2] << endl;
				}
				//U[t + 1].S[0][ne] -= uu[0];
				//U[t + 1].S[1][ne] -= uu[1] / uu[0];
				//U[t + 1].S[2][ne] -= uu[2] / uu[0];

				//U[t + 1].S[0][nk] -= normal.x * (U[t].S[1][nk] + U[t].S[2][ne]) / 2 + normal.y * (U[t].S[2][nk] + U[t].S[2][ne]) / 2 +
				//lambda_x / 2 * (U[t].S[0][ne] - U[t].S[0][nk]);
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
				//vector<double> wk(3);
				//wk[0] = U[t].S[0][nk];
				//wk[1] = U[t].S[1][nk] * U[t].S[0][nk];
				//wk[2] = U[t].S[2][nk] * U[t].S[0][nk];

				//lambda_x = U[t].S[1][nk] + sqrt(g * U[t].S[0][nk]);
				//lambda_y = U[t].S[2][nk] + sqrt(g * U[t].S[0][nk]);

				//lambdax = std::max(lambdax, lambda_x);
				//lambday = std::max(lambday, lambda_y);

				//vector<double> uk(3);
				//uk[0] = U[t].S[0][nk];
				//uk[1] = U[t].S[1][nk];
				//uk[2] = U[t].S[2][nk];

				//const auto rk = R(uk);
				//const auto gk = G(uk);

				//W[t + 1].S[0][nk] -= t_step / dx * (normal.x * (rk[0]) + normal.y * (gk[0]) + lambda_x / 2 * (uk[0] - uk[0]));
				//W[t + 1].S[1][nk] -= t_step / dx * (normal.x * (rk[1]) + normal.y * (gk[1]) + lambda_x / 2 * (-uk[1] - uk[1]));
				//W[t + 1].S[2][nk] -= t_step / dx * (normal.x * (rk[2]) + normal.y * (gk[2]) + lambda_x / 2 * (-uk[2] - uk[2]));

				//U[t + 1].S[0][nk] -= W[t + 1].S[0][nk];
				//U[t + 1].S[1][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];
				//U[t + 1].S[2][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];

				vector<double> u(3);
				u[0] = U[t].S[0][nk];
				u[1] = -U[t].S[1][nk];
				u[2] = -U[t].S[2][nk];

				//lambda_x = U[t].S[1][nk] / U[t].S[1][nk] + sqrt(g * U[t].S[0][nk]);
				//lambda_y = U[t].S[2][nk] / U[t].S[0][nk] + sqrt(g * U[t].S[0][ne]);

				//lambda_x = std::max(fabs(U[t].S[1][nk]) + sqrt(g * U[t].S[0][nk]), fabs(u[1]) + sqrt(g * u[0]));
				//lambda_y = std::max(fabs(U[t].S[2][nk]) + sqrt(g * U[t].S[0][nk]), fabs(u[2]) + sqrt(g * u[0]));

				lambdax = std::max(lambdax, lambda_x);
				lambday = std::max(lambday, lambda_y);

				vector<double> uk(3);
				uk[0] = U[t].S[0][nk];
				uk[1] = U[t].S[1][nk];
				uk[2] = U[t].S[2][nk];

				vector<double> ue(3);
				ue[0] = u[0];
				ue[1] = u[1];
				ue[2] = u[2];

				const auto rk = R(uk);
				const auto re = R(ue);
				const auto gk = G(uk);
				const auto ge = G(ue);
				lambda_x = std::max(fabs(uk[1]) + sqrt(g * uk[0]), fabs(ue[1]) + sqrt(g * ue[0]));
				lambda_y = std::max(fabs(uk[2]) + sqrt(g * uk[0]), fabs(ue[2]) + sqrt(g * ue[0]));
				vector<double> uu(3);
				//if (nk == 4)
				//	normal.x = -normal.x;
				if (normal.x > 0 || normal.y > 0)
				{
					uu[0] = t_step / dx * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
					uu[1] = t_step / dx * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
					uu[2] = t_step / dx * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
				}
				else
				{
					uu[0] = t_step / dx * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
					uu[1] = t_step / dx * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
					uu[2] = t_step / dx * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
				}



				if (nk == 7)
				{
					/*cout << "nn" << endl;
					cout << uu[1] << endl;
					cout << lambda_x << endl;
					cout << lambda_y << endl;
					cout << ue[1] << endl;
					cout << uk[1] << endl;
					cout << normal.x << endl;
					cout << normal.y << endl;
					cout << normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 << endl;
					cout << (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])) << endl;*/
					//cout << uu[1] << endl;
					//cout << uu[2] << endl;
				}

				W[t + 1].S[0][nk] -= uu[0];
				W[t + 1].S[1][nk] -= uu[1];
				W[t + 1].S[2][nk] -= uu[2];

				//U[t + 1].S[0][nk] -= W[t + 1].S[0][nk];
				//U[t + 1].S[1][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];
				//U[t + 1].S[2][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];

				//W[t + 1].S[0][ne] -= -normal.x * (rk[0] + re[0]) / 2 - normal.y * (gk[0] + ge[0]) / 2 + lambda_x / 2 * (ue[0] - uk[0]);
				//W[t + 1].S[1][ne] -= -normal.x * (rk[1] + re[1]) / 2 - normal.y * (gk[1] + ge[1]) / 2 + lambda_x / 2 * (ue[1] - uk[1]);
				//W[t + 1].S[2][ne] -= -normal.x * (rk[2] + re[2]) / 2 - normal.y * (gk[2] + ge[2]) / 2 + lambda_x / 2 * (ue[2] - uk[2]);

				//U[t + 1].S[0][ne] -= W[t + 1].S[0][nk];
				//U[t + 1].S[1][ne] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];
				//U[t + 1].S[2][ne] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];


				//U[t + 1].S[0][nk] = U[t].S[0][nk];
				//U[t + 1].S[1][nk] = -U[t].S[1][nk];
				//U[t + 1].S[2][nk] = -U[t].S[2][nk];

				//W[t + 1].S[0][nk] = W[t].S[0][nk];
				//W[t + 1].S[1][nk] = -W[t].S[1][nk];
				//W[t + 1].S[2][nk] = -W[t].S[2][nk];
			}
		}
		for (size_t i = 0; i < size; ++i)
		{

			U[t + 1].S[0][i] = W[t + 1].S[0][i];
			U[t + 1].S[1][i] = W[t + 1].S[1][i] / W[t + 1].S[0][i];
			U[t + 1].S[2][i] = W[t + 1].S[2][i] / W[t + 1].S[0][i];
		}
		W.push_back(solution(size));
		U.push_back(solution(size));
		
	}
	return 0;
}

const int dg_solver_shallow_water::solve(
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
	const function<const vector<double>(const vector<double>&)>& R,
	const function<const vector<double>(const vector<double>&)>& G,
	const function<const vector<double>(const vector<double>&)>& F) const
{
	// U = [H Hu Hv]
	// R = [Hu Hu^2 Huv]
	// G = [Hv Huv Hv^2]
	// S = [0 -gHdkdx+1/p(t_sx-t_bx+Fx -gHdkdy+1/p(t_sy-t_by+Fy))]
	// test
	// U = [h hu hv]
	// R = [hu hu^2+gh^2/2 huv]
	// G = [hv huv hv^2+gh^2/2]
	// S = [0 0 0]
	vector<double> Ut[3];
	const int max_iter = 30000;
	const double dx = (x1 - x0) / nx;
	const double dy = (y1 - y0) / ny;
    CRegularMesh mesh{ x0, y0, x1, y1, (int)nx, (int)ny };
	// dunno
	auto center = [=](const int i)
	{
		const auto& elem = mesh.GetElement(i);
		vector<Point> pts(4);
		pts[0] = mesh.GetNode(elem->GetNode(0));
		pts[1] = mesh.GetNode(elem->GetNode(1));
		pts[2] = mesh.GetNode(elem->GetNode(2));
		pts[3] = mesh.GetNode(elem->GetNode(3));
		return Point(pts[0].x + (pts[3].x - pts[0].x) / 2, pts[0].y + (pts[3].y - pts[0].y) / 2);
	};
	const int size = mesh.GetNumberOfElements();
	const int bsize = mesh.GetNumberOfBoundaries();
	struct solution
	{
		vector<double> S[3];
		solution() {}
		solution(const int _size)
		{
			S[0].resize(_size);
			S[1].resize(_size);
			S[2].resize(_size);
		}
	};
	vector<solution> U(2);
	vector<solution> W(2);
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
		const auto cent(center(i));
		if (cent.x >= -1. / 2 && cent.x <= 1. / 2 && cent.y >= -1. / 2 && cent.y <= 1. / 2)
		{
			U[0].S[0][i] = 4;
			U[0].S[1][i] = 0;
			U[0].S[2][i] = 0;

			W[0].S[0][i] = 4;
			W[0].S[1][i] = 0;
			W[0].S[2][i] = 0;
		}
		else
		{
			U[0].S[0][i] = 1;
			U[0].S[1][i] = 0;
			U[0].S[2][i] = 0;

			W[0].S[0][i] = 1;
			W[0].S[1][i] = 0;
			W[0].S[2][i] = 0;
		}
	}
	//parameter<double> g(1);
	
	//parameter<double> H(0);
	
	double t_step = 0.1;
	const double cfl = 0.8;
	// W = [h hu hv]
	double lambda_x = 0;
	double lambda_y = 0;
	double lambdax = 0;
	double lambday = 0;
	double t_curr = 0;
	size_t iter_max = 10000;
	for (size_t t = 0; t < iter_max && t_curr < t1; ++t, t_curr += t_step)
	{
		cout << "dt:\t" << t_step << endl;
		cout << "t:\t" << t_curr << endl;
		lambda_x = 0;
		lambda_y = 0;
		for (size_t i = 0; i < size; ++i)
		{
			const auto& res = F(vector<double>{W[t].S[0][i], W[t].S[1][i], W[t].S[2][i]});
			W[t + 1].S[0][i] = W[t].S[0][i] + res[0];
			W[t + 1].S[1][i] = W[t].S[1][i] + res[1];// *U[t].S[1][i];
			W[t + 1].S[2][i] = W[t].S[2][i] + res[2];// *U[t].S[2][i];

											//U[t + 1].S[0][i] = U[t].S[0][i];
											//U[t + 1].S[1][i] = U[t].S[1][i];
											//U[t + 1].S[2][i] = U[t].S[2][i];

			lambda_x = std::max(fabs(U[t].S[1][i]) + sqrt(g * U[t].S[0][i]), lambda_x);
			lambda_y = std::max(fabs(U[t].S[2][i]) + sqrt(g * U[t].S[0][i]), lambda_y);
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
				vector<double> wk(3);
				wk[0] = U[t].S[0][nk];
				wk[1] = U[t].S[1][nk] * U[t].S[0][nk];
				wk[2] = U[t].S[2][nk] * U[t].S[0][nk];

				vector<double> we(3);
				we[0] = U[t].S[0][ne];
				we[1] = U[t].S[1][ne] * U[t].S[0][ne];
				we[2] = U[t].S[2][ne] * U[t].S[0][ne];

				//lambda_x = U[t].S[1][nk] / U[t].S[1][nk] + sqrt(g * U[t].S[0][nk]);
				//lambda_y = U[t].S[2][nk] / U[t].S[0][nk] + sqrt(g * U[t].S[0][ne]);

				lambda_x = std::max(fabs(U[t].S[1][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[1][ne]) + sqrt(g * U[t].S[0][ne]));
				lambda_y = std::max(fabs(U[t].S[2][nk]) + sqrt(g * U[t].S[0][nk]), fabs(U[t].S[2][ne]) + sqrt(g * U[t].S[0][ne]));

				lambdax = std::max(lambdax, lambda_x);
				lambday = std::max(lambday, lambda_y);
				double ll = max(lambda_x, lambda_y);
				//cout << "max:\t" << ll << endl;
				vector<double> uk(3);
				uk[0] = U[t].S[0][nk];
				uk[1] = U[t].S[1][nk];
				uk[2] = U[t].S[2][nk];

				vector<double> ue(3);
				ue[0] = U[t].S[0][ne];
				ue[1] = U[t].S[1][ne];
				ue[2] = U[t].S[2][ne];

				const auto rk = R(uk);
				const auto re = R(ue);
				const auto gk = G(uk);
				const auto ge = G(ue);

				vector<double> uu(3);
				uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[0] + re[0]) / 2 + normal.y * (gk[0] + ge[0]) / 2 - (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
				uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 - (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
				uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (normal.x * (rk[2] + re[2]) / 2 + normal.y * (gk[2] + ge[2]) / 2 - (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
				if (nk == 7)
				{
					//cout << "nk" << endl;
					//cout << uu[1] << endl;
					//cout << uu[1] << endl;
					//cout << uu[2] << endl;
				}
				W[t + 1].S[0][nk] -= uu[0];
				W[t + 1].S[1][nk] -= uu[1];
				W[t + 1].S[2][nk] -= uu[2];

				//U[t + 1].S[0][nk] -= uu[0];
				//U[t + 1].S[1][nk] -= uu[1] / uu[0];
				//U[t + 1].S[2][nk] -= uu[2] / uu[0];

				uu[0] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[0] + re[0]) / 2 - normal.y * (gk[0] + ge[0]) / 2 + (lambda_x * normal.x / 2 * (ue[0] - uk[0]) + lambda_y * normal.y / 2 * (ue[0] - uk[0])));
				uu[1] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[1] + re[1]) / 2 - normal.y * (gk[1] + ge[1]) / 2 + (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])));
				uu[2] = t_step / mesh.GetElement(nk)->GetMeasure() * bound->GetMeasure() * (-normal.x * (rk[2] + re[2]) / 2 - normal.y * (gk[2] + ge[2]) / 2 + (lambda_x * normal.x / 2 * (ue[2] - uk[2]) + lambda_y * normal.y / 2 * (ue[2] - uk[2])));
				W[t + 1].S[0][ne] -= uu[0];
				W[t + 1].S[1][ne] -= uu[1];
				W[t + 1].S[2][ne] -= uu[2];
				if (ne == 7)
				{
					//cout << "ne" << endl;
					//cout << uu[1] << endl;
					//cout << uu[1] << endl;
					//cout << uu[2] << endl;
				}
				//U[t + 1].S[0][ne] -= uu[0];
				//U[t + 1].S[1][ne] -= uu[1] / uu[0];
				//U[t + 1].S[2][ne] -= uu[2] / uu[0];

				//U[t + 1].S[0][nk] -= normal.x * (U[t].S[1][nk] + U[t].S[2][ne]) / 2 + normal.y * (U[t].S[2][nk] + U[t].S[2][ne]) / 2 +
				//lambda_x / 2 * (U[t].S[0][ne] - U[t].S[0][nk]);
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
				//vector<double> wk(3);
				//wk[0] = U[t].S[0][nk];
				//wk[1] = U[t].S[1][nk] * U[t].S[0][nk];
				//wk[2] = U[t].S[2][nk] * U[t].S[0][nk];

				//lambda_x = U[t].S[1][nk] + sqrt(g * U[t].S[0][nk]);
				//lambda_y = U[t].S[2][nk] + sqrt(g * U[t].S[0][nk]);

				//lambdax = std::max(lambdax, lambda_x);
				//lambday = std::max(lambday, lambda_y);

				//vector<double> uk(3);
				//uk[0] = U[t].S[0][nk];
				//uk[1] = U[t].S[1][nk];
				//uk[2] = U[t].S[2][nk];

				//const auto rk = R(uk);
				//const auto gk = G(uk);

				//W[t + 1].S[0][nk] -= t_step / dx * (normal.x * (rk[0]) + normal.y * (gk[0]) + lambda_x / 2 * (uk[0] - uk[0]));
				//W[t + 1].S[1][nk] -= t_step / dx * (normal.x * (rk[1]) + normal.y * (gk[1]) + lambda_x / 2 * (-uk[1] - uk[1]));
				//W[t + 1].S[2][nk] -= t_step / dx * (normal.x * (rk[2]) + normal.y * (gk[2]) + lambda_x / 2 * (-uk[2] - uk[2]));

				//U[t + 1].S[0][nk] -= W[t + 1].S[0][nk];
				//U[t + 1].S[1][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];
				//U[t + 1].S[2][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];

				vector<double> u(3);
				u[0] = U[t].S[0][nk];
				u[1] = -U[t].S[1][nk];
				u[2] = -U[t].S[2][nk];

				//lambda_x = U[t].S[1][nk] / U[t].S[1][nk] + sqrt(g * U[t].S[0][nk]);
				//lambda_y = U[t].S[2][nk] / U[t].S[0][nk] + sqrt(g * U[t].S[0][ne]);

				//lambda_x = std::max(fabs(U[t].S[1][nk]) + sqrt(g * U[t].S[0][nk]), fabs(u[1]) + sqrt(g * u[0]));
				//lambda_y = std::max(fabs(U[t].S[2][nk]) + sqrt(g * U[t].S[0][nk]), fabs(u[2]) + sqrt(g * u[0]));

				lambdax = std::max(lambdax, lambda_x);
				lambday = std::max(lambday, lambda_y);

				vector<double> uk(3);
				uk[0] = U[t].S[0][nk];
				uk[1] = U[t].S[1][nk];
				uk[2] = U[t].S[2][nk];

				vector<double> ue(3);
				ue[0] = u[0];
				ue[1] = u[1];
				ue[2] = u[2];

				const auto rk = R(uk);
				const auto re = R(ue);
				const auto gk = G(uk);
				const auto ge = G(ue);
				lambda_x = std::max(fabs(uk[1]) + sqrt(g * uk[0]), fabs(ue[1]) + sqrt(g * ue[0]));
				lambda_y = std::max(fabs(uk[2]) + sqrt(g * uk[0]), fabs(ue[2]) + sqrt(g * ue[0]));
				vector<double> uu(3);
				//if (nk == 4)
				//	normal.x = -normal.x;
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



				if (nk == 7)
				{
					/*cout << "nn" << endl;
					cout << uu[1] << endl;
					cout << lambda_x << endl;
					cout << lambda_y << endl;
					cout << ue[1] << endl;
					cout << uk[1] << endl;
					cout << normal.x << endl;
					cout << normal.y << endl;
					cout << normal.x * (rk[1] + re[1]) / 2 + normal.y * (gk[1] + ge[1]) / 2 << endl;
					cout << (lambda_x * normal.x / 2 * (ue[1] - uk[1]) + lambda_y * normal.y / 2 * (ue[1] - uk[1])) << endl;*/
					//cout << uu[1] << endl;
					//cout << uu[2] << endl;
				}

				W[t + 1].S[0][nk] -= uu[0];
				W[t + 1].S[1][nk] -= uu[1];
				W[t + 1].S[2][nk] -= uu[2];

				//U[t + 1].S[0][nk] -= W[t + 1].S[0][nk];
				//U[t + 1].S[1][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];
				//U[t + 1].S[2][nk] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];

				//W[t + 1].S[0][ne] -= -normal.x * (rk[0] + re[0]) / 2 - normal.y * (gk[0] + ge[0]) / 2 + lambda_x / 2 * (ue[0] - uk[0]);
				//W[t + 1].S[1][ne] -= -normal.x * (rk[1] + re[1]) / 2 - normal.y * (gk[1] + ge[1]) / 2 + lambda_x / 2 * (ue[1] - uk[1]);
				//W[t + 1].S[2][ne] -= -normal.x * (rk[2] + re[2]) / 2 - normal.y * (gk[2] + ge[2]) / 2 + lambda_x / 2 * (ue[2] - uk[2]);

				//U[t + 1].S[0][ne] -= W[t + 1].S[0][nk];
				//U[t + 1].S[1][ne] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];
				//U[t + 1].S[2][ne] -= W[t + 1].S[0][nk] / W[t + 1].S[0][nk];


				//U[t + 1].S[0][nk] = U[t].S[0][nk];
				//U[t + 1].S[1][nk] = -U[t].S[1][nk];
				//U[t + 1].S[2][nk] = -U[t].S[2][nk];

				//W[t + 1].S[0][nk] = W[t].S[0][nk];
				//W[t + 1].S[1][nk] = -W[t].S[1][nk];
				//W[t + 1].S[2][nk] = -W[t].S[2][nk];
			}
		}
		for (size_t i = 0; i < size; ++i)
		{

			U[t + 1].S[0][i] = W[t + 1].S[0][i];
			U[t + 1].S[1][i] = W[t + 1].S[1][i] / W[t + 1].S[0][i];
			U[t + 1].S[2][i] = W[t + 1].S[2][i] / W[t + 1].S[0][i];
		}
		W.push_back(solution(size));
		U.push_back(solution(size));

	}
	return 0;
}




