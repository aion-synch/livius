#include "test_case_elliptic_fem.h"
#include "../CoreNCFEM/Grids/TriangularMesh.h"
#include "../CoreNCFEM/Grids/RegularMesh.h"
#include "../CoreNCFEM/Methods/FEMethod.h"
#include "../Problems/DiffusionScalar.h"
#include "../CoreNCA/MatrixSkyline.h"
#include "../CoreNCFEM/Methods/FEAnalysis.h"
#include "../Solvers/fem_solver.h"
#include "../CoreNCFEM/GaussianField.h"
#include "../CoreNCFEM/FiniteElements/Triangle.h"
#include <random>
#define _USE_MATH_DEFINES

#include <math.h>
using namespace corenc;
using namespace std;
using namespace Mesh;
using namespace Algebra;
using namespace method;
test_case_elliptic_fem::test_case_elliptic_fem()
{

}

test_case_elliptic_fem::~test_case_elliptic_fem()
{

}

const double kekus(const double c, const double a = 0, const double b = 90.)
{
	return a + b + c;
}

const int test_case_elliptic_fem::elliptic_fem_square_lin_basis() const
{
	const auto x0 = [=](const Point& p)
	{
		return 1.;
	};
	const auto x1 = [=](const Point& p)
	{
		//return 10. - p.x - p.y;
		return (10. - p.x) * (10. - p.y);
	};
	const auto x11 = [=](const Point& p)
	{
		return 10. - p.x - p.y - p.x * p.y;
	};
	const auto x2 = [=](const Point& p)
	{
		//return 10. - p.x - p.y - p.x * p.x - p.y * p.y;
		return (10. - p.x - p.y) * (10. - p.x - p.y);
	};
	const auto x3 = [=](const Point& p)
	{
		return (10. - p.x - p.y) * (10. - p.x - p.y) * (10. - p.x - p.y);
	};
	const auto x4 = [=](const Point& p)
	{
		return (10. - p.x - p.y) * (10. - p.x - p.y) * (10. - p.x - p.y) * (10. - p.x - p.y);
		//return std::pow(10. - p.x - p.y, 40);
		//return exp((1. - p.x) * (1. - p.y));
	};
	const auto x5 = [=](const Point& p)
	{
		return (10. - p.x - p.y) * (10. - p.x - p.y) * (10. - p.x - p.y) * (10. - p.x - p.y) * (10. - p.x - p.y);
	};
	const auto xexp = [=](const Point& p)
	{
		return exp((1. - p.x) * (1. - p.y));
	};
	const auto f1 = [=](const Point& p)
	{
		return 0.;
	};
	const auto f2 = [=](const Point& p)
	{
		return -4.;
	};
	const auto f3 = [=](const Point& p)
	{
		return 12. * (p.x + p.y) - 120.;
	};
	const auto f4 = [=](const Point& p)
	{
		return -24. * (p.x + p.y - 10.) * (p.x + p.y - 10.);
		//return -3120. * std::pow((p.x + p.y - 10.), 38);
		//return -exp((p.x - 1)*(p.y - 1))*(p.x - 1) *(p.x - 1) - exp((p.x - 1)*(p.y - 1))*(p.y - 1) * (p.y - 1);
	};
	const auto f5 = [=](const Point& p)
	{
		return -24. * (p.x + p.y - 10.) * (p.x + p.y - 10.);
	};
	const auto fexp = [=](const Point& p)
	{
		return -exp((p.x - 1) * (p.y - 1)) * ((p.x - 1) * (p.x - 1) + (p.y - 1) * (p.y - 1));
	};
	
	
	CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
	mesh.refine_p();
	mesh.refine_p();
	//mesh.refine_p();
	const auto x1b = [=](const int el, const int node, const Point& p, const std::function<const double(const Point&)>& func)
	{
		if(node < 3)
			return func(p);
		const auto& elem = mesh.GetElement(el);
		vector<Point> pts{ mesh.GetNode(elem->GetNode(0)), mesh.GetNode(elem->GetNode(1)), mesh.GetNode(elem->GetNode(2)) };
		const double q3 = 4 * func(Point{ pts[0].x + (pts[1].x - pts[0].x) / 2, pts[0].y + (pts[1].y - pts[0].y) / 2 }) - 2 * (func(pts[0]) + func(pts[1]));
		const double q4 = 4 * func(Point{ pts[1].x + (pts[2].x - pts[1].x) / 2, pts[1].y + (pts[2].y - pts[1].y) / 2 }) - 2 * (func(pts[1]) + func(pts[2]));
		const double q5 = 4 * func(Point{ pts[0].x + (pts[2].x - pts[0].x) / 2, pts[0].y + (pts[2].y - pts[0].y) / 2 }) - 2 * (func(pts[0]) + func(pts[2]));
		const Point center = Point{ (pts[0].x + pts[1].x + pts[2].x) / 3, (pts[0].y + pts[1].y + pts[2].y) / 3 };
		const double q6 = 27. / 2 * func(Point{ pts[0].x + (pts[1].x - pts[0].x) / 3, pts[0].y + (pts[1].y - pts[0].y) / 3 }) - 9 * func(pts[0]) - 9. / 2 * func(pts[1]) - 3 * q3;
		const double q7 = 27. / 2 * func(Point{ pts[1].x + (pts[2].x - pts[1].x) / 3, pts[1].y + (pts[2].y - pts[1].y) / 3 }) - 9 * func(pts[1]) - 9. / 2 * func(pts[2]) - 3 * q4;
		const double q8 = 27. / 2 * func(Point{ pts[0].x + (pts[2].x - pts[0].x) / 3, pts[0].y + (pts[2].y - pts[0].y) / 3 }) - 9 * func(pts[0]) - 9. / 2 * func(pts[2]) - 3 * q5;
		const double q9 = 27. * (func(center) - 1. / 3 * (func(pts[0]) + func(pts[1]) + func(pts[2])) - 1. / 9 * (q3 + q4 + q5));

		/*const Point tp(0.15, 0.5);
		const double expected = func(tp);
		const double actual = GetShapeFunction(0, tp) * func(pts[0])
		+ GetShapeFunction(1, tp) * func(pts[1])
		+ GetShapeFunction(2, tp) * func(pts[2])
		+ GetShapeFunction(3, tp) * q3
		+ GetShapeFunction(4, tp) * q4
		+ GetShapeFunction(5, tp) * q5
		+ GetShapeFunction(6, tp) * q6
		+ GetShapeFunction(7, tp) * q7
		+ GetShapeFunction(8, tp) * q8
		+ GetShapeFunction(9, tp) * q9;*/
		switch (node)
		{
		case 3:
			return q3;
		case 4:
			return q4;
		case 5:
			return q5;
		case 6:
			return q6;
		}
		return 0.;
	};
	const auto bdn = [=](const int el, const int node, const Point& p)
	{
		//return x1b(el, node, p, x2);
		return x1(p);
	};
	const auto src = [=](const int el, const int node, const Point& p)
	{
		//return 4.;
		return f1(p);
		//return x1b(el, node, p, f3);
	};
	const parameter<double> boundary_lin(bdn);
	const parameter<double> source(src);

	CDiffusionScalar problem;
	problem.addTerm(Terms::EFV);
	problem.add_parameter(Terms::EFV, 6, source);
	problem.add_parameter(Terms::IDUDV, 6, 1.);
	problem.add_boundary_parameter(1, 1, boundary_lin);
	problem.add_boundary_parameter(1, 2, boundary_lin);
	problem.add_boundary_parameter(1, 3, boundary_lin);
	problem.add_boundary_parameter(1, 4, boundary_lin);
	solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>> fem;
	vector<double> solution;
	vector<double> exact(mesh.GetNumberOfNodes());
	vector<double> exact2(mesh.GetNumberOfNodes());
	vector<double> actual(exact);
	fem.elliptic_solver(&problem, &mesh, &solution);
	for (int i = 0; i < exact.size(); ++i)
	{
		exact[i] = x4(mesh.GetNode(i));
		actual[i] = fem.get_value(mesh, solution, mesh.GetNode(i));
	}
	const auto& N_el = mesh.GetNumberOfElements();
	double sum = 0.;
	double sum2 = 0.;
	for (int i = 0; i < N_el; ++i)
	{
		const auto& elem = mesh.GetElement(i);
		vector<Point> pts{ mesh.GetNode(elem->GetNode(0)), mesh.GetNode(elem->GetNode(1)), mesh.GetNode(elem->GetNode(2)) };
		//elem->GetWeight(8, pts, x4);
		//exact2[elem->GetNode(0)] = x1b(i, 0, pts[0], x2);
		//exact2[elem->GetNode(1)] = x1b(i, 1, pts[1], x2);
		//exact2[elem->GetNode(2)] = x1b(i, 2, pts[2], x2);
		//exact2[elem->GetNode(3)] = x1b(i, 3, mesh.GetNode(elem->GetNode(3)), x2);
		//exact2[elem->GetNode(4)] = x1b(i, 4, mesh.GetNode(elem->GetNode(4)), x2);
		//exact2[elem->GetNode(5)] = x1b(i, 5, mesh.GetNode(elem->GetNode(5)), x2);
		sum += elem->Integrate([&](const Point& p){
			const double temp = fem.get_value(mesh, solution, p);
			return (temp - x1(p)) * (temp - x1(p)); }, pts);
		sum2 += elem->Integrate([&](const Point& p){return x1(p) * x1(p); }, pts);
	}
	ofstream ofs1("log.txt");
	const double sum1 = sqrt(sum);
	ofs1.precision(7);
	scientific(ofs1);
	ofs1 << "Absolute L2-norm: " << sum1 << endl;
	sum = sqrt(sum / sum2);
	ofs1 << "Relative L2-norm: " << sum << endl;
	return 0;
}

const int test_case_elliptic_fem::elliptic_fem_hp_fixed(const int h_ref_max, const int p_ref_max) const
{
	ofstream ofs1("log.txt");
	std::random_device rd;
	std::mt19937 mt(rd());
	//std::uniform_real_distribution<double> distribution(, -4);
	const double a_max = 1;
	const int N = 10;
	std::normal_distribution<double> distribution(0, a_max);
	std::uniform_real_distribution<double> distribution_uni(0, 1);
	std::vector<double> sum(2);
	vector<double> temp(N);
	vector<Point> centers(N);
	ofs1 << "Generated Random Values" << endl;
	for (auto i = 0; i < N; ++i)
	{
		temp[i] = distribution(mt);
		centers[i].x = distribution_uni(mt);
		centers[i].y = distribution_uni(mt);
		ofs1 << temp[i] << endl;
	}
	ofs1 << "Centers X" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].x << endl;
	ofs1 << "Centers Y" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].y << endl;
	GaussianKernel gk{ N, centers };
	for (auto i = 0; i < 10; ++i)
		for (auto j = 0; j < 10; ++j)
		{
			cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
		}
	cout.precision(15);
	scientific(cout);
	ofs1.precision(15);
	scientific(ofs1);
	CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
	vector<vector<CTriangularMesh>> meshes;
	vector<vector<vector<double>>> solutions;
	vector<vector<double>> norms;
	vector<vector<double>> abs_norms;
	vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
	meshes.resize(p_ref_max);
	for (auto& it : meshes)
		it.resize(h_ref_max);
	fems.resize(p_ref_max);
	solutions.resize(p_ref_max);
	norms.resize(p_ref_max);
	abs_norms.resize(p_ref_max);
	for (auto i = 0; i < p_ref_max; ++i)
	{
		fems[i].resize(h_ref_max);
		solutions[i].resize(h_ref_max);
		norms[i].resize(h_ref_max);
		abs_norms[i].resize(h_ref_max);
		for (auto j = 0; j < h_ref_max; ++j)
		{
			meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
			for (auto ii = 0; ii < i; ++ii)
				meshes[i][j].refine_p();
		}
	}
	const auto fexp = [=](const Point& p)
	{
		return 100.;
		return 100 * exp(-10 * ((p.x - 0.5) * (p.x - 0.5) + (p.y - 0.5) * (p.y - 0.5)));
		//return std::pow((p.x + p.y - 10.), 4);
		//return -3120. * std::pow((p.x + p.y - 10.), 38);
		//return -exp((p.x - 1) * (p.y - 1)) * ((p.x - 1) * (p.x - 1) + (p.y - 1) * (p.y - 1));
		//return -4.*(p.x - 10)*(p.y - 10) - (p.x - 10)*(2.*p.x + 2.*p.y - 20.) - (p.y - 10.)*(2.*p.x + 2.*p.y - 20.);
		//return -10.*std::pow(p.x+p.y-10.,4);
		//return -2.*exp(10. - p.x - p.y);
		//return 80.*exp(10 - p.x - p.y)*std::pow(p.x + p.y - 10., 39) - 3120.*exp(10. - p.y - p.x)*std::pow(p.x + p.y - 10., 38);
		//return exp(1 - p.x - p.y)*exp((p.x - 1) * (p.y - 1)) * (p.x - 1) + exp(1. - p.x - p.y) * exp((p.x - 1.) * (p.y - 1)) * (p.y - 1) - exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.x - 1) * (p.x - 1) -
		//	exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.y - 1) * (p.y - 1);
	};
	const auto bdn = [=](const int el, const int node, const Point& p)
	{
		//return exp((1 - p.x)*(1 - p.y));
		//return std::pow(10. - p.x - p.y, 40);
		return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto bdn1 = [=](const int el, const int node, const Point& p)
	{
		//return exp((1 - p.x)*(1 - p.y));
		//return std::pow(10. - p.x - p.y, 40);
		return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto src = [=](const int el, const int node, const Point& p)
	{
		return fexp(p);
	};

	double _max = 0;
	double _min = 1000;
	const auto lam = [&](const int el, const int node, const Point& p)
	{
		//return exp(1. - p.x - p.y);
		//return std::pow(10. - p.x - p.y, 40);
		//return 1.;
		const double val = exp(gk.get_gp(temp, p));
		//const double val = gk.get_gp(temp, p);
		if (_max < val)
		{
			_max = val;
			cout << "max: " << val << endl;
		}
		if (_min > val)
		{
			_min = val;
			cout << "min: " << val << endl;
		}
		return val;
	};
	const parameter<double> boundary_lin(bdn);
	const parameter<double> boundary_lin1(bdn1);
	const parameter<double> source(src);
	const parameter<double> lamda(lam);
	CDiffusionScalar problem;
	problem.addTerm(Terms::EFV);
	problem.add_parameter(Terms::EFV, 6, source);
	problem.add_parameter(Terms::IDUDV, 6, lamda);
	problem.add_boundary_parameter(1, 1, boundary_lin);
	problem.add_boundary_parameter(1, 2, boundary_lin);
	problem.add_boundary_parameter(1, 3, boundary_lin1);
	problem.add_boundary_parameter(1, 4, boundary_lin);
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				norms[i][j] += elem->Integrate([&](const Point& p) 
				{
					const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p); return t * t;
				}, pts);
			}
			norms[i][j] = sqrt(norms[i][j]);
			ofs1 << i << j << ": \t" << norms[i][j] << endl;
		}
	}
	const int Nx = 100;
	const int Ny = 100;
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			ofstream res("results//mesh" + to_string(i) + to_string(j) + ".txt");
			for (auto ii = 0; ii < Nx + 1; ++ii)
			{
				for (auto jj = 0; jj < Ny + 1; ++jj)
				{
					const auto xx = (double)ii / Nx;
					const auto yy = (double)jj / Ny;
					res << xx << "\t" << yy << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii/Nx, (double)jj/Ny)) << endl;
				}
			}
			//const auto& N_el = meshes[i][j].GetNumberOfElements();
		}
	}
	ofs1 << "Relative L2-norm" << endl;
	for (auto i = 0; i < p_ref_max - 1; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
						fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
						//bdn(ii, 0, p);
					return t * t; 
				}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
		}
	}

	for (auto i = p_ref_max - 1; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max - 1; ++j)
		{
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
						fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
					//bdn(ii, 0, p);
					return t * t;
				}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
		}
	}
	ofs1 << "Convergence order" << endl;
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max - 1; ++j)
		{
			ofs1 << i << j << ": \t" << std::log2(abs_norms[i][j] / abs_norms[i][j + 1]) << endl;
		}
	}
	return 0;
}

const int test_case_elliptic_fem::elliptic_fem_hp_fixed_triangle(const int h_ref_max, const int p_ref_max) const
{
	ofstream ofs1("log.txt");
	std::random_device rd;
	std::mt19937 mt(rd());
	//std::uniform_real_distribution<double> distribution(, -4);
	const double a_max = 1;
	const int N = 10;
	std::normal_distribution<double> distribution(0, a_max);
	std::uniform_real_distribution<double> distribution_uni(0, 1);
	std::vector<double> sum(2);
	vector<double> temp(N);
	vector<Point> centers(N);
	ofs1 << "Generated Random Values" << endl;
	for (auto i = 0; i < N; ++i)
	{
		temp[i] = distribution(mt);
		centers[i].x = distribution_uni(mt);
		centers[i].y = distribution_uni(mt);
		ofs1 << temp[i] << endl;
	}
	ofs1 << "Centers X" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].x << endl;
	ofs1 << "Centers Y" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].y << endl;
	GaussianKernel gk{ N, centers };
	for (auto i = 0; i < 10; ++i)
		for (auto j = 0; j < 10; ++j)
		{
			cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
		}
	cout.precision(15);
	scientific(cout);
	ofs1.precision(15);
	scientific(ofs1);
	//CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
	CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2};
    vector<vector<CTriangularMesh>> meshes;
	vector<vector<vector<double>>> solutions;
	vector<vector<double>> norms;
	vector<vector<double>> abs_norms;
	vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
	meshes.resize(p_ref_max);
	for (auto& it : meshes)
		it.resize(h_ref_max);
	fems.resize(p_ref_max);
	solutions.resize(p_ref_max);
	norms.resize(p_ref_max);
	abs_norms.resize(p_ref_max);
	for (auto i = 0; i < p_ref_max; ++i)
	{
		fems[i].resize(h_ref_max);
		solutions[i].resize(h_ref_max);
		norms[i].resize(h_ref_max);
		abs_norms[i].resize(h_ref_max);
		for (size_t j = 0; j < h_ref_max; ++j)
		{
			//meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
            meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j), (int)std::pow(2, j)};
			for (auto ii = 0; ii < i; ++ii)
                meshes[i][j].refine_p();
		}
	}
	const auto fexp = [=](const Point& p)
	{
        return 100. * (p.x * p.x + p.y * p.y) * cos(10. * p.x * p.y);
        //return -exp(p.x * p.y) * (p.x * p.x * p.x * p. y + 2 * p.x * p.x + p.x * p.y * p.y * p.y + 2 * p.y * p.y);
        //return 0.;
        return -12. * (p.x * p.x + p.y * p.y);
		//return 100.;
		//return 100 * exp(-10 * ((p.x - 0.5) * (p.x - 0.5) + (p.y - 0.5) * (p.y - 0.5)));
		//return std::pow((p.x + p.y - 10.), 4);
		return -3120. * std::pow((p.x + p.y - 10.), 38);
		//return -exp((p.x - 1) * (p.y - 1)) * ((p.x - 1) * (p.x - 1) + (p.y - 1) * (p.y - 1));
		//return -4.*(p.x - 10)*(p.y - 10) - (p.x - 10)*(2.*p.x + 2.*p.y - 20.) - (p.y - 10.)*(2.*p.x + 2.*p.y - 20.);
		//return -10.*std::pow(p.x+p.y-10.,4);
		//return -2.*exp(10. - p.x - p.y);
		//return 80.*exp(10 - p.x - p.y)*std::pow(p.x + p.y - 10., 39) - 3120.*exp(10. - p.y - p.x)*std::pow(p.x + p.y - 10., 38);
		//return exp(1 - p.x - p.y)*exp((p.x - 1) * (p.y - 1)) * (p.x - 1) + exp(1. - p.x - p.y) * exp((p.x - 1.) * (p.y - 1)) * (p.y - 1) - exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.x - 1) * (p.x - 1) -
		//	exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.y - 1) * (p.y - 1);
	};
	const auto bdn = [=](const int el, const int node, const Point& p)
	{
        return cos(10. * p.x * p.y);
        //return p.x * p.y * exp(p.x * p.y);
        //return 10.;
        return p.x * p.x * p.x * p.x + p.y * p.y * p.y * p.y;
		//return exp((1 - p.x)*(1 - p.y));
		return std::pow(10. - p.x - p.y, 40);
		//return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto bdn1 = [=](const int el, const int node, const Point& p)
	{
        return cos(10. * p.x * p.y);
        //return p.x * p.y * exp(p.x * p.y);
        //return 10.;
        return p.x * p.x * p.x * p.x + p.y * p.y * p.y * p.y;
		//return exp((1 - p.x)*(1 - p.y));
		return std::pow(10. - p.x - p.y, 40);
		//return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto src = [=](const int el, const int node, const Point& p)
	{
		return fexp(p);
	};

	double _max = 0;
	double _min = 1000;
	const auto lam = [&](const int el, const int node, const Point& p)
	{
		//return exp(1. - p.x - p.y);
		//return std::pow(10. - p.x - p.y, 40);
		return 1.;
		const double val = exp(gk.get_gp(temp, p));
		//const double val = gk.get_gp(temp, p);
		if (_max < val)
		{
			_max = val;
			cout << "max: " << val << endl;
		}
		if (_min > val)
		{
			_min = val;
			cout << "min: " << val << endl;
		}
		return val;
	};
	const parameter<double> boundary_lin(bdn);
	const parameter<double> boundary_lin1(bdn1);
	const parameter<double> source(src);
	const parameter<double> lamda(lam);
	CDiffusionScalar problem;
	problem.addTerm(Terms::EFV);
	problem.add_parameter(Terms::EFV, 0, source);
	problem.add_parameter(Terms::IDUDV, 0, lamda);
	problem.add_boundary_parameter(1, 0, boundary_lin);
	problem.add_boundary_parameter(1, 1, boundary_lin);
	problem.add_boundary_parameter(1, 2, boundary_lin1);
	problem.add_boundary_parameter(1, 3, boundary_lin);
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
            fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
            const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				norms[i][j] += elem->Integrate([&](const Point& p)
					{
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p); return t * t;
					}, pts);
			}
			norms[i][j] = sqrt(norms[i][j]);
			ofs1 << i << j << ": \t" << norms[i][j] << endl;
		}
	}
	const int Nx = 100;
	const int Ny = 100;
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			ofstream res("results//mesh" + to_string(i) + to_string(j) + ".txt");
			for (auto ii = 0; ii < Nx + 1; ++ii)
			{
				for (auto jj = 0; jj < Ny + 1; ++jj)
				{
					const auto xx = (double)ii / Nx;
					const auto yy = (double)jj / Ny;
                    res << xx << "\t" << yy << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii / Nx, (double)jj / Ny)) << endl;
				}
			}
			//const auto& N_el = meshes[i][j].GetNumberOfElements();
		}
	}
	ofs1 << "Relative L2-norm" << endl;
	for (auto i = 0; i < p_ref_max - 1; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
            const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
                        fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						//fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
					bdn(ii, 0, p);
					return t * t;
					}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
		}
	}

	for (auto i = p_ref_max - 1; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max - 1; ++j)
		{
            const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
                        fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						//fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
					bdn(ii, 0, p);
					return t * t;
					}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
		}
	}
	ofs1 << "Convergence order for h" << endl;
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max - 1; ++j)
		{
            ofs1 << i << j << ": \t" << std::log2(abs_norms[i][j] / abs_norms[i][j + 1]) << endl;
		}
	}
	ofs1 << "Convergence order for p" << endl;
	for (auto i = 0; i < h_ref_max; ++i)
	{
		for (auto j = 0; j < p_ref_max - 1; ++j)
		{
			ofs1 << j << i << ": \t" << std::log2(abs_norms[j][i] / abs_norms[j + 1][i]) << endl;
		}
    }
	return 0;
}


const int test_case_elliptic_fem::elliptic_fem_hp_lagrange_triangle(const int h_ref_max, const int p_ref_max) const
{
    ofstream ofs1("log.txt");
    std::random_device rd;
    std::mt19937 mt(rd());
    //std::uniform_real_distribution<double> distribution(, -4);
    const double a_max = 1;
    const int N = 10;
    std::normal_distribution<double> distribution(0, a_max);
    std::uniform_real_distribution<double> distribution_uni(0, 1);
    std::vector<double> sum(2);
    vector<double> temp(N);
    vector<Point> centers(N);
    ofs1 << "Generated Random Values" << endl;
    for (auto i = 0; i < N; ++i)
    {
        temp[i] = distribution(mt);
        centers[i].x = distribution_uni(mt);
        centers[i].y = distribution_uni(mt);
        ofs1 << temp[i] << endl;
    }
    ofs1 << "Centers X" << endl;
    for (auto i = 0; i < N; ++i)
        ofs1 << centers[i].x << endl;
    ofs1 << "Centers Y" << endl;
    for (auto i = 0; i < N; ++i)
        ofs1 << centers[i].y << endl;
    GaussianKernel gk{ N, centers };
    for (auto i = 0; i < 10; ++i)
        for (auto j = 0; j < 10; ++j)
        {
            cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
        }
    cout.precision(15);
    scientific(cout);
    ofs1.precision(15);
    scientific(ofs1);
    //CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
    CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2};
    vector<vector<CTriangularMesh>> meshes;
    vector<vector<vector<double>>> solutions;
    vector<vector<double>> norms;
    vector<vector<double>> abs_norms;
    vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
    meshes.resize(p_ref_max);
    for (auto& it : meshes)
        it.resize(h_ref_max);
    fems.resize(p_ref_max);
    solutions.resize(p_ref_max);
    norms.resize(p_ref_max);
    abs_norms.resize(p_ref_max);
    for (auto i = 0; i < p_ref_max; ++i)
    {
        fems[i].resize(h_ref_max);
        solutions[i].resize(h_ref_max);
        norms[i].resize(h_ref_max);
        abs_norms[i].resize(h_ref_max);
        for (size_t j = 0; j < h_ref_max; ++j)
        {
            //meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
            meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j), (int)std::pow(2, j)};
            //for (auto ii = 0; ii < i; ++ii)
                meshes[i][j].set4thOrder();
        }
    }
    cout << "test" << endl;
    const auto fexp = [=](const Point& p)
    {
        return 100 * (p.x * p.x + p.y * p.y) * cos(10. * p.x * p.y);
        //return p.x + p.y;
        //return 0.;
        //return -12 * (p.x * p.x + p.y * p.y);
        return -exp(p.x * p.y) * (p.x * p.x * p.x * p. y + 2 * p.x * p.x + p.x * p.y * p.y * p.y + 2 * p.y * p.y);
        //return -3120. * std::pow((p.x + p.y - 10.), 38);
    };
    const auto bdn = [=](const int el, const int node, const Point& p)
    {
        return cos(10 * p.x * p.y);
        //return 10.;
        //return p.x + p.y;
        //return p.x * p.x * p.x * p.x + p.y * p.y * p.y * p.y;
        return p.x * p.y * exp(p.x * p.y);
        return std::pow(10. - p.x - p.y, 40);
    };
    const auto bdn1 = [=](const int el, const int node, const Point& p)
    {
        return cos(10 * p.x * p.y);
        //return p.x * p.x * p.x * p.x + p.y * p.y * p.y * p.y;
        //return 10.;
        //return p.x + p.y;
        return p.x * p.y * exp(p.x * p.y);
        return std::pow(10. - p.x - p.y, 40);
    };
    const auto src = [=](const int el, const int node, const Point& p)
    {
        return fexp(p);
    };

    double _max = 0;
    double _min = 1000;
    const auto lam = [&](const int el, const int node, const Point& p)
    {
        return 1.;
        const double val = exp(gk.get_gp(temp, p));
        if (_max < val)
        {
            _max = val;
            cout << "max: " << val << endl;
        }
        if (_min > val)
        {
            _min = val;
            cout << "min: " << val << endl;
        }
        return val;
    };
    const parameter<double> boundary_lin(bdn);
    const parameter<double> boundary_lin1(bdn);
    const parameter<double> source(src);
    const parameter<double> lamda(lam);
    CDiffusionScalar problem;
    problem.addTerm(Terms::EFV);
    //problem.addTerm(Terms::IUV);
    problem.add_parameter(Terms::EFV, 0, source);
    problem.add_parameter(Terms::IDUDV, 0, lamda);
    //problem.add_parameter(Terms::IUV, 0, lamda);
    problem.add_boundary_parameter(1, 0, boundary_lin);
    problem.add_boundary_parameter(1, 1, boundary_lin);
    problem.add_boundary_parameter(1, 2, boundary_lin1);
    problem.add_boundary_parameter(1, 3, boundary_lin);
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max; ++j)
        {
            fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
            //for (int ii = 0; ii < meshes[i][j].GetNumberOfNodes(); ++ii)
            //{
              //  cout << bdn(0, 0, meshes[i][j].GetNode(ii)) << "\t" << solutions[i][j][ii] << endl;
            //}
            const auto& N_el = meshes[i][j].GetNumberOfElements();
            for (auto ii = 0; ii < N_el; ++ii)
            {
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                norms[i][j] += elem->Integrate([&](const Point& p)
                    {
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p); return t * t;
                    }, pts);
            }
            norms[i][j] = sqrt(norms[i][j]);
            ofs1 << i << j << ": \t" << norms[i][j] << endl;
        }
    }
    const int Nx = 100;
    const int Ny = 100;
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max; ++j)
        {
            ofstream res("results//mesh" + to_string(i) + to_string(j) + ".txt");
            for (auto ii = 0; ii < Nx + 1; ++ii)
            {
                for (auto jj = 0; jj < Ny + 1; ++jj)
                {
                    const auto xx = (double)ii / Nx;
                    const auto yy = (double)jj / Ny;
                    res << xx << "\t" << yy << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii / Nx, (double)jj / Ny)) << endl;
                }
            }
            //const auto& N_el = meshes[i][j].GetNumberOfElements();
        }
    }
    ofs1 << "Relative L2-norm" << endl;
    for (auto i = 0; i < p_ref_max - 1; ++i)
    {
        for (auto j = 0; j < h_ref_max; ++j)
        {
            const auto& N_el = meshes[i][j].GetNumberOfElements();
            for (auto ii = 0; ii < N_el; ++ii)
            {
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                abs_norms[i][j] += elem->Integrate([&](const Point& p) {
                    const double t =
                        fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
                        //fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
                    bdn(ii, 0, p);
                    return t * t;
                    }, pts);
            }
            abs_norms[i][j] = sqrt(abs_norms[i][j]);
            ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
            ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
        }
    }

    for (auto i = p_ref_max - 1; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max - 1; ++j)
        {
            const auto& N_el = meshes[i][j].GetNumberOfElements();
            for (auto ii = 0; ii < N_el; ++ii)
            {
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                abs_norms[i][j] += elem->Integrate([&](const Point& p) {
                    const double t =
                        fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
                        //fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
                    bdn(ii, 0, p);
                    return t * t;
                    }, pts);
            }
            abs_norms[i][j] = sqrt(abs_norms[i][j]);
            ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
            ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
        }
    }
    ofs1 << "Convergence order for h" << endl;
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max - 1; ++j)
        {
            ofs1 << i << j << ": \t" << std::log2(abs_norms[i][j] / abs_norms[i][j + 1]) << endl;
        }
    }
    ofs1 << "Convergence order for p" << endl;
    for (auto i = 0; i < h_ref_max; ++i)
    {
        for (auto j = 0; j < p_ref_max - 1; ++j)
        {
            ofs1 << j << i << ": \t" << std::log2(abs_norms[j][i] / abs_norms[j + 1][i]) << endl;
        }
    }
    return 0;
}

const int test_case_elliptic_fem::elliptic_fem_hxhy_fixed_triangle(const int hx_max, const int hy_max) const
{
	ofstream ofs1("log.txt");
	std::random_device rd;
	std::mt19937 mt(rd());
	//std::uniform_real_distribution<double> distribution(, -4);
	const double a_max = 1;
	const int N = 10;
	std::normal_distribution<double> distribution(0, a_max);
	std::uniform_real_distribution<double> distribution_uni(0, 1);
	std::vector<double> sum(2);
	vector<double> temp(N);
	vector<Point> centers(N);
	ofs1 << "Generated Random Values" << endl;
	for (auto i = 0; i < N; ++i)
	{
		temp[i] = distribution(mt);
		centers[i].x = distribution_uni(mt);
		centers[i].y = distribution_uni(mt);
		ofs1 << temp[i] << endl;
	}
	ofs1 << "Centers X" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].x << endl;
	ofs1 << "Centers Y" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].y << endl;
	GaussianKernel gk{ N, centers };
	for (auto i = 0; i < 10; ++i)
		for (auto j = 0; j < 10; ++j)
		{
			cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
		}
	cout.precision(15);
	scientific(cout);
	ofs1.precision(15);
	scientific(ofs1);
	//CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
	CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2 };
	vector<vector<CTriangularMesh>> meshes;
	vector<vector<vector<double>>> solutions;
	vector<vector<double>> norms;
	vector<vector<double>> abs_norms;
	vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
	meshes.resize(hy_max);
	for (auto& it : meshes)
		it.resize(hx_max);
	fems.resize(hy_max);
	solutions.resize(hy_max);
	norms.resize(hy_max);
	abs_norms.resize(hy_max);
	for (auto i = 0; i < hy_max; ++i)
	{
		fems[i].resize(hx_max);
		solutions[i].resize(hx_max);
		norms[i].resize(hx_max);
		abs_norms[i].resize(hx_max);
		for (size_t j = 0; j < hx_max; ++j)
		{
			//meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
			meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j + 1), (int)std::pow(2, i + 1) };
		}
	}
	const auto fexp = [=](const Point& p)
	{
		//return 100.;
		return 100 * exp(-10 * ((p.x - 0.5) * (p.x - 0.5) + (p.y - 0.5) * (p.y - 0.5)));
		//return std::pow((p.x + p.y - 10.), 4);
		return -3120. * std::pow((p.x + p.y - 10.), 38);
		//return -exp((p.x - 1) * (p.y - 1)) * ((p.x - 1) * (p.x - 1) + (p.y - 1) * (p.y - 1));
		//return -4.*(p.x - 10)*(p.y - 10) - (p.x - 10)*(2.*p.x + 2.*p.y - 20.) - (p.y - 10.)*(2.*p.x + 2.*p.y - 20.);
		//return -10.*std::pow(p.x+p.y-10.,4);
		//return -2.*exp(10. - p.x - p.y);
		//return 80.*exp(10 - p.x - p.y)*std::pow(p.x + p.y - 10., 39) - 3120.*exp(10. - p.y - p.x)*std::pow(p.x + p.y - 10., 38);
		//return exp(1 - p.x - p.y)*exp((p.x - 1) * (p.y - 1)) * (p.x - 1) + exp(1. - p.x - p.y) * exp((p.x - 1.) * (p.y - 1)) * (p.y - 1) - exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.x - 1) * (p.x - 1) -
		//	exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.y - 1) * (p.y - 1);
	};
	const auto bdn = [=](const int el, const int node, const Point& p)
	{
		//return exp((1 - p.x)*(1 - p.y));
		//return std::pow(10. - p.x - p.y, 40);
		return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto bdn1 = [=](const int el, const int node, const Point& p)
	{
		//return exp((1 - p.x)*(1 - p.y));
		//return std::pow(10. - p.x - p.y, 40);
		return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto src = [=](const int el, const int node, const Point& p)
	{
		return fexp(p);
	};

	double _max = 0;
	double _min = 1000;
	const auto lam = [&](const int el, const int node, const Point& p)
	{
		//return exp(1. - p.x - p.y);
		//return std::pow(10. - p.x - p.y, 40);
		//return 1.;
		const double val = exp(gk.get_gp(temp, p));
		//const double val = gk.get_gp(temp, p);
		if (_max < val)
		{
			_max = val;
			cout << "max: " << val << endl;
		}
		if (_min > val)
		{
			_min = val;
			cout << "min: " << val << endl;
		}
		return val;
	};
	const parameter<double> boundary_lin(bdn);
	const parameter<double> boundary_lin1(bdn1);
	const parameter<double> source(src);
	const parameter<double> lamda(lam);
	CDiffusionScalar problem;
	problem.addTerm(Terms::EFV);
	problem.add_parameter(Terms::EFV, 0, source);
	problem.add_parameter(Terms::IDUDV, 0, lamda);
	problem.add_boundary_parameter(1, 0, boundary_lin);
	problem.add_boundary_parameter(1, 1, boundary_lin);
	problem.add_boundary_parameter(1, 2, boundary_lin1);
	problem.add_boundary_parameter(1, 3, boundary_lin);
	for (auto i = 0; i < hy_max; ++i)
	{
		for (auto j = 0; j < hx_max; ++j)
		{
			fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				norms[i][j] += elem->Integrate([&](const Point& p)
					{
						const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p); return t * t;
					}, pts);
			}
			norms[i][j] = sqrt(norms[i][j]);
			ofs1 << i << j << ": \t" << norms[i][j] << endl;
		}
	}
	const int Nx = 100;
	const int Ny = 100;
	for (auto i = 0; i < hy_max; ++i)
	{
		for (auto j = 0; j < hx_max; ++j)
		{
			ofstream res("results//mesh" + to_string(i) + to_string(j) + ".txt");
			for (auto ii = 0; ii < Nx + 1; ++ii)
			{
				for (auto jj = 0; jj < Ny + 1; ++jj)
				{
					const auto xx = (double)ii / Nx;
					const auto yy = (double)jj / Ny;
					res << xx << "\t" << yy << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii / Nx, (double)jj / Ny)) << endl;
				}
			}
			//const auto& N_el = meshes[i][j].GetNumberOfElements();
		}
	}
	ofs1 << "Relative L2-norm" << endl;
	for (auto i = 0; i < hy_max - 1; ++i)
	{
		for (auto j = 0; j < hx_max; ++j)
		{
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
						fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						fems[hy_max - 1][hx_max - 1].get_value(meshes[hy_max - 1][hx_max - 1], solutions[hy_max - 1][hx_max - 1], p);
						//bdn(ii, 0, p);
					return t * t;
					}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[hy_max - 1][hx_max - 1] << endl;
		}
	}

	for (auto i = hy_max - 1; i < hy_max; ++i)
	{
		for (auto j = 0; j < hx_max - 1; ++j)
		{
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
						fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						fems[hy_max - 1][hx_max - 1].get_value(meshes[hy_max - 1][hx_max - 1], solutions[hy_max - 1][hx_max - 1], p);
						//bdn(ii, 0, p);
					return t * t;
					}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[hy_max - 1][hx_max - 1] << endl;
		}
	}
	ofs1 << "Convergence order for hx" << endl;
	for (auto i = 0; i < hy_max; ++i)
	{
		for (auto j = 0; j < hx_max - 1; ++j)
		{
			ofs1 << i << j << ": \t" << std::log2(abs_norms[i][j] / abs_norms[i][j + 1]) << endl;
		}
	}
	ofs1 << "Convergence order for hy" << endl;
	for (auto i = 0; i < hx_max; ++i)
	{
		for (auto j = 0; j < hy_max - 1; ++j)
		{
			ofs1 << j << i << ": \t" << std::log2(abs_norms[j][i] / abs_norms[j + 1][i]) << endl;
		}
	}
	return 0;
}


const int test_case_elliptic_fem::global_matrix(const int h_ref_max, const int p_ref_max) const
{
    ofstream ofs1("log.txt");
    std::random_device rd;
    std::mt19937 mt(rd());
    //std::uniform_real_distribution<double> distribution(, -4);
    const double a_max = 1;
    const int N = 10;
    std::normal_distribution<double> distribution(0, a_max);
    std::uniform_real_distribution<double> distribution_uni(0, 1);
    std::vector<double> sum(2);
    vector<double> temp(N);
    vector<Point> centers(N);
    ofs1 << "Generated Random Values" << endl;
    for (auto i = 0; i < N; ++i)
    {
        temp[i] = distribution(mt);
        centers[i].x = distribution_uni(mt);
        centers[i].y = distribution_uni(mt);
        ofs1 << temp[i] << endl;
    }
    ofs1 << "Centers X" << endl;
    for (auto i = 0; i < N; ++i)
        ofs1 << centers[i].x << endl;
    ofs1 << "Centers Y" << endl;
    for (auto i = 0; i < N; ++i)
        ofs1 << centers[i].y << endl;
    GaussianKernel gk{ N, centers };
    for (auto i = 0; i < 10; ++i)
        for (auto j = 0; j < 10; ++j)
        {
            cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
        }
    cout.precision(15);
    scientific(cout);
    ofs1.precision(15);
    scientific(ofs1);
    //CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
    CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2};
    vector<vector<CTriangularMesh>> meshes;
    vector<vector<vector<double>>> solutions;
    vector<vector<double>> norms;
    vector<vector<double>> abs_norms;
    vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
    meshes.resize(p_ref_max);
    for (auto& it : meshes)
        it.resize(h_ref_max);
    fems.resize(p_ref_max);
    solutions.resize(p_ref_max);
    norms.resize(p_ref_max);
    abs_norms.resize(p_ref_max);
    for (auto i = 0; i < p_ref_max; ++i)
    {
        fems[i].resize(h_ref_max);
        solutions[i].resize(h_ref_max);
        norms[i].resize(h_ref_max);
        abs_norms[i].resize(h_ref_max);
        for (size_t j = 0; j < h_ref_max; ++j)
        {
            //meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
            meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j), (int)std::pow(2, j)};
            //for (auto ii = 0; ii < i; ++ii)
                meshes[i][j].set4thOrder();
        }
    }
    const auto fexp = [=](const Point& p)
    {
        //return p.x + p.y;
        //return 10.;
        return 0.;
        return -12 * (p.x * p.x + p.y * p.y);
        return -exp(p.x * p.y) * (p.x * p.x * p.x * p. y + 2 * p.x * p.x + p.x * p.y * p.y * p.y + 2 * p.y * p.y);
        //return -3120. * std::pow((p.x + p.y - 10.), 38);
    };
    const auto bdn = [=](const int el, const int node, const Point& p)
    {
        //return 10.;
        return p.x + p.y;
        return p.x * p.x * p.x * p.x + p.y * p.y * p.y * p.y;
        return p.x * p.y * exp(p.x * p.y);
        return std::pow(10. - p.x - p.y, 40);
    };
    const auto bdn1 = [=](const int el, const int node, const Point& p)
    {
        //return 10.;
        return p.x + p.y;
        return p.x * p.y * exp(p.x * p.y);
        return std::pow(10. - p.x - p.y, 40);
    };
    const auto src = [=](const int el, const int node, const Point& p)
    {
        return fexp(p);
    };

    double _max = 0;
    double _min = 1000;
    const auto lam = [&](const int el, const int node, const Point& p)
    {
        return 1.;
        const double val = exp(gk.get_gp(temp, p));
        if (_max < val)
        {
            _max = val;
            cout << "max: " << val << endl;
        }
        if (_min > val)
        {
            _min = val;
            cout << "min: " << val << endl;
        }
        return val;
    };
    const parameter<double> boundary_lin(bdn);
    const parameter<double> boundary_lin1(bdn);
    const parameter<double> source(src);
    const parameter<double> lamda(lam);
    CDiffusionScalar problem;
    problem.addTerm(Terms::EFV);
    //problem.addTerm(Terms::IUV);
    problem.add_parameter(Terms::EFV, 0, source);
    problem.add_parameter(Terms::IDUDV, 0, lamda);
    //problem.add_parameter(Terms::IUV, 0, lamda);
    problem.add_boundary_parameter(1, 0, boundary_lin);
    problem.add_boundary_parameter(1, 1, boundary_lin);
    problem.add_boundary_parameter(1, 2, boundary_lin1);
    problem.add_boundary_parameter(1, 3, boundary_lin);
    for (int k = 0; k < h_ref_max; ++k)
    //vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
    {
        cout << "Number of elements: " << meshes[0][k].GetNumberOfElements() << endl;
        for (int i = 0; i < meshes[0][k].GetNumberOfElements(); ++i)
        {
            auto elem = meshes[0][k].GetElement(i);
            for (int j = 0; j < elem->GetDoFs(); ++j)
                cout << elem->GetNode(j) << "\t";
            cout << endl;
        }
        cout << "Number of nodes: " << meshes[0][k].GetNumberOfNodes() << endl;
        for (int i = 0; i < meshes[0][k].GetNumberOfNodes(); ++i)
        {
            cout << i << "\t" << meshes[0][k].GetNode(i).x << "\t" << meshes[0][k].GetNode(i).y << endl;
        }
        std::vector<double> res;
        std::vector<double> res2;
        //std::shared_ptr<Algebra::MatrixSkyline> matrix{ new Algebra::MatrixSkyline() };
        //Algebra::MatrixSkyline* matrix{ new Algebra::MatrixSkyline() };
        Algebra::Matrix* matrix{ new Algebra::Matrix() };
        std::vector<double> rhs;
        auto m_method = new FEMethod<CDiffusionScalar, CTriangularMesh, Algebra::Matrix>{ &problem, &meshes[0][k], matrix, &rhs };
        cout << "Discretization" << endl;
        m_method->Discretization();
        std::vector<double> tts(matrix->GetSize());
        Algebra::ESolver esl{ Algebra::Solvers::BiCGStab };
        //solutions[0][h_ref_max - 1] = esl.Solve(*matrix, rhs, solutions[0][h_ref_max - 1], res, 100000, 1e-13);
        solutions[0][k].resize(matrix->GetSize());
        cout << "Solving" << endl;
        esl.Gauss(*matrix, rhs, solutions[0][k]);
        cout << "Checking" << endl;
        for (int i = 0; i < matrix->GetSize(); ++i)
        {
            for (int j = 0; j < matrix->GetSize(); ++j)
            {
                tts[i] += matrix->operator()(i, j) * bdn(0, 0, meshes[0][k].GetNode(j));
                cout << matrix->operator()(i, j) << "\t";
            }

            cout << rhs[i] << endl;
            if (fabs(tts[i] - rhs[i]) > 1e-12)
                cout << "Something wrong with me: " << fabs(tts[i] - rhs[i]) << endl;
            cout << solutions[0][k][i] << "\t" << bdn(0, 0, meshes[0][k].GetNode(i)) << endl;
        }
        delete m_method;
        delete matrix;
    }
    return 0;
}


const int test_case_elliptic_fem::conv_diff_fem_fixed_triangle(const int h_ref_max, const int p_ref_max) const
{
	ofstream ofs1("log.txt");
	std::random_device rd;
	std::mt19937 mt(rd());
	//std::uniform_real_distribution<double> distribution(, -4);
	const double a_max = 1;
	const int N = 10;
	std::normal_distribution<double> distribution(0, a_max);
	std::uniform_real_distribution<double> distribution_uni(0, 1);
	std::vector<double> sum(2);
	vector<double> temp(N);
	vector<Point> centers(N);
	ofs1 << "Generated Random Values" << endl;
	for (auto i = 0; i < N; ++i)
	{
		temp[i] = distribution(mt);
		centers[i].x = distribution_uni(mt);
		centers[i].y = distribution_uni(mt);
		ofs1 << temp[i] << endl;
	}
	ofs1 << "Centers X" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].x << endl;
	ofs1 << "Centers Y" << endl;
	for (auto i = 0; i < N; ++i)
		ofs1 << centers[i].y << endl;
	GaussianKernel gk{ N, centers };
	for (auto i = 0; i < 10; ++i)
		for (auto j = 0; j < 10; ++j)
		{
			cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
		}
	cout.precision(15);
	scientific(cout);
	ofs1.precision(15);
	scientific(ofs1);
	//CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
	CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2 };
	vector<vector<CTriangularMesh>> meshes;
	vector<vector<vector<double>>> solutions;
	vector<vector<double>> norms;
	vector<vector<double>> abs_norms;
	vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
	meshes.resize(p_ref_max);
	for (auto& it : meshes)
		it.resize(h_ref_max);
	fems.resize(p_ref_max);
	solutions.resize(p_ref_max);
	norms.resize(p_ref_max);
	abs_norms.resize(p_ref_max);
    for (auto i = 0; i < p_ref_max; ++i)
	{
		fems[i].resize(h_ref_max);
		solutions[i].resize(h_ref_max);
		norms[i].resize(h_ref_max);
		abs_norms[i].resize(h_ref_max);
		for (size_t j = 0; j < h_ref_max; ++j)
		{
			//meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
            meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j + 1), (int)std::pow(2, j + 1) };
			for (auto ii = 0; ii < i; ++ii)
				meshes[i][j].refine_p();
		}
	}
    const double lal = 1e0;
    const double v = 1e6;
	const auto fexp = [=](const Point& p)
    {
        return v * (p.y * exp(p.x * p.y) + p.x * exp(p.x * p.y)) - lal * (p.y * p.y + p.x * p.x) * exp(p.x * p.y);
        return 2. * v;
		//return v;
		//return 20 * p.x + 2 * p.y + 0.4;
		//return 100.;
		//return 100 * exp(-10 * ((p.x - 0.5) * (p.x - 0.5) + (p.y - 0.5) * (p.y - 0.5)));
		//return std::pow((p.x + p.y - 10.), 4);
        //return -3120. * std::pow((p.x + p.y - 10.), 38);
		//return -exp((p.x - 1) * (p.y - 1)) * ((p.x - 1) * (p.x - 1) + (p.y - 1) * (p.y - 1));
		//return -4.*(p.x - 10)*(p.y - 10) - (p.x - 10)*(2.*p.x + 2.*p.y - 20.) - (p.y - 10.)*(2.*p.x + 2.*p.y - 20.);
		//return -10.*std::pow(p.x+p.y-10.,4);
		//return -2.*exp(10. - p.x - p.y);
		//return 80.*exp(10 - p.x - p.y)*std::pow(p.x + p.y - 10., 39) - 3120.*exp(10. - p.y - p.x)*std::pow(p.x + p.y - 10., 38);
		//return exp(1 - p.x - p.y)*exp((p.x - 1) * (p.y - 1)) * (p.x - 1) + exp(1. - p.x - p.y) * exp((p.x - 1.) * (p.y - 1)) * (p.y - 1) - exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.x - 1) * (p.x - 1) -
		//	exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.y - 1) * (p.y - 1);
	};
	const auto bdn = [=](const int el, const int node, const Point& p)
	{
        return exp(p.x * p.y);
		return p.x + p.y;
		return p.x * p.x + p.y * p.y;
		//return exp((1 - p.x)*(1 - p.y));
		//return std::pow(10. - p.x - p.y, 40);
		//return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto bdn1 = [=](const int el, const int node, const Point& p)
	{
		//return exp((1 - p.x)*(1 - p.y));
		return std::pow(10. - p.x - p.y, 40);
		//return 0.;
		//return std::pow(10. - p.x - p.y, 2);
		//return 10. - p.x - p.y;
	};
	const auto src = [=](const int el, const int node, const Point& p)
	{
		return fexp(p);
	};

	double _max = 0;
	double _min = 1000;
	const auto lam = [&](const int el, const int node, const Point& p)
	{
		//return exp(1. - p.x - p.y);
		//return std::pow(10. - p.x - p.y, 40);
        return lal;
		const double val = exp(gk.get_gp(temp, p));
		//const double val = gk.get_gp(temp, p);
		if (_max < val)
		{
			_max = val;
			cout << "max: " << val << endl;
		}
		if (_min > val)
		{
			_min = val;
			cout << "min: " << val << endl;
		}
		return val;
	};
	const auto vel = [&](const int el, const int node, const Point& p)
	{
		return Point(v, v);
	};
	const parameter<double> boundary_lin(bdn);
	const parameter<double> boundary_lin1(bdn1);
	const parameter<double> source(src);
	const parameter<double> lamda(lam);
	const parameter<Point>	velocity(vel);
	CDiffusionScalar problem;
	problem.addTerm(Terms::EFV);
    problem.add_parameter(Terms::EFV, 0, source);
    problem.addTerm(Terms::IDUV);
    problem.add_parameter(Terms::IDUV, 0, velocity);
    problem.add_parameter(Terms::IDUDV, 0, lamda);
	problem.add_boundary_parameter(1, 0, boundary_lin);
	problem.add_boundary_parameter(1, 1, boundary_lin);
	problem.add_boundary_parameter(1, 2, boundary_lin);
	problem.add_boundary_parameter(1, 3, boundary_lin);
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				norms[i][j] += elem->Integrate([&](const Point& p)
					{
						const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p); return t * t;
					}, pts);
			}
			norms[i][j] = sqrt(norms[i][j]);
			ofs1 << i << j << ": \t" << norms[i][j] << endl;
		}
	}
	const int Nx = 100;
	const int Ny = 100;
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			ofstream res("results//mesh" + to_string(i) + to_string(j) + ".txt");
			for (auto ii = 0; ii < Nx + 1; ++ii)
			{
				for (auto jj = 0; jj < Ny + 1; ++jj)
				{
					const auto xx = (double)ii / Nx;
					const auto yy = (double)jj / Ny;
					res << xx << "\t" << yy << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii / Nx, (double)jj / Ny)) << endl;
				}
			}
			//const auto& N_el = meshes[i][j].GetNumberOfElements();
		}
	}
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max; ++j)
        {
            ofstream res("results//line" + to_string(i) + to_string(j) + ".txt");
            for (auto ii = 0; ii < Nx + 1; ++ii)
            {
                //for (auto jj = 0; jj < Ny + 1; ++jj)
                {
                    const auto xx = (double)ii / Nx;
                    //const auto yy = (double)jj / Ny;
                    res << xx << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii / Nx, 0.5)) << endl;
                }
            }
            //const auto& N_el = meshes[i][j].GetNumberOfElements();
        }
    }
    /*ofs1 << "Relative L2-norm" << endl;
    for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max; ++j)
		{
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
						fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						//fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
						bdn(ii, 0, p);
					return t * t;
					}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
		}
	}

	for (auto i = p_ref_max - 1; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max - 1; ++j)
		{
			const auto& N_el = meshes[i][j].GetNumberOfElements();
			for (auto ii = 0; ii < N_el; ++ii)
			{
				const auto& elem = meshes[i][j].GetElement(ii);
				vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
				abs_norms[i][j] += elem->Integrate([&](const Point& p) {
					const double t =
						fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
						//fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
						bdn(ii, 0, p);
					return t * t;
					}, pts);
			}
			abs_norms[i][j] = sqrt(abs_norms[i][j]);
			ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
			ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
		}
	}
	ofs1 << "Convergence order for h" << endl;
	for (auto i = 0; i < p_ref_max; ++i)
	{
		for (auto j = 0; j < h_ref_max - 1; ++j)
		{
			ofs1 << i << j << ": \t" << std::log2(abs_norms[i][j] / abs_norms[i][j + 1]) << endl;
		}
	}
	ofs1 << "Convergence order for p" << endl;
	for (auto i = 0; i < h_ref_max; ++i)
	{
		for (auto j = 0; j < p_ref_max - 1; ++j)
		{
			ofs1 << j << i << ": \t" << std::log2(abs_norms[j][i] / abs_norms[j + 1][i]) << endl;
		}
    }*/
	return 0;
}




const int test_case_elliptic_fem::elliptic_fem_2d_tria() const
{
	solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>> fem;
	cout << "Performing the test \"2d elliptic via fem on triangles\"" << endl;
	CTriangularMesh mesh2{ "grids//2d//unit_square_layer//unit.msh" };
	shared_ptr<CDiffusionScalar> problem{ new CDiffusionScalar() };
	vector<double> res;
	vector<double> dg_solution;
	vector<vector<double>> all_souls;
	vector<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>> meths;
	vector<CTriangularMesh> meshes;
	vector<double> conv;
	if (problem->load_parameters("grids//2d//unit_square_layer//param.txt"))
		return 1;
	{
		shared_ptr<MatrixSkyline> matrix{ new MatrixSkyline() };
		vector<double> rhs; 
		CTriangularMesh mesh{ "grids//2d//unit_square_layer//unit.msh" };
		vector<double> solution{ mesh.getSolution() };
		FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline> method{ problem.get(), &mesh, matrix.get(), &rhs };
		method.Discretization();
		ESolver esl{ Solvers::GMRES };
		solution = esl.Solve(*matrix.get(), rhs, solution, res, 100000, 1e-13);
		all_souls.push_back(solution);
		meths.push_back(method);
		meshes.push_back(mesh);
		fem.elliptic_solver(problem.get(), &mesh, &solution);
	}
	
	cout << "-----------------------------" << endl;
	cout << "h-refinement test convergence" << endl;
	cout << "-----------------------------" << endl;
	conv.push_back(0.);

	for (int i = 0; i < 4; ++i)
	{
		vector<double> rhs2;
		shared_ptr<MatrixSkyline> matrix2{ new MatrixSkyline() };
		mesh2.refine_h();
		cout << "refine:\t" << i << endl;
		cout << "Number of the nodes:\t" << mesh2.GetNumberOfNodes() << endl;
		vector<double> solution2{ mesh2.getSolution() };
		FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline> method2{ problem.get(), &mesh2, matrix2.get(), &rhs2 };
		method2.Discretization();
		ESolver esl2{ Solvers::GMRES };
		solution2 = esl2.Solve(*matrix2.get(), rhs2, solution2, res, 100000, 1e-13);
		meths.push_back(method2);
		meshes.push_back(mesh2);
		FEAnalysis<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, CTriangularMesh, CTriangularMesh> anal;

		conv.push_back(anal.L2Norm(meths[i], method2, meshes[i], mesh2, all_souls[i], solution2));
		cout << "Relative convergence in L2-norm:\t" << conv[i + 1] << endl;
		//cout << "Absolute convergence in L2-norm:\t" << anal.L2Norm(meths[0], method2, meshes[0], mesh2, all_souls[0], solution2) << endl;
		cout << "Convergence rate:\t" << (fabs(conv[i] / conv[i + 1])) << "\t not log" << endl;
		cout << "Convergence rate:\t" << log2(fabs(conv[i] / conv[i + 1])) << endl;
		all_souls.push_back(solution2);
		cout << method2.GetSolution(mesh2, solution2, Point(0.4, 0.4, 0.)) << endl;
	}
	double t2 = 0;
	cout << "-----------------------------" << endl;
	for (int i = 0; i < 4; ++i)
	{
		FEAnalysis<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, CTriangularMesh, CTriangularMesh> anal;
		const auto temp = anal.L2Norm(meths[i], meths[3], meshes[i], meshes[3], all_souls[i], all_souls[3]);
		cout << "Absolute convergence in L2-norm:\t" << temp << endl;
		cout << "Convergence rate:\t" << log2(fabs(t2 / temp)) << endl;
		t2 = temp;
	}
	CTriangularMesh mesh3{ "grids//2d//unit_square_layer//unit.msh" };

	cout << "-----------------------------" << endl;
	cout << "p-refinement test convergence" << endl;
	cout << "-----------------------------" << endl;

	for (int i = 0; i < 6; ++i)
	{
		ofstream ifs("p" + to_string(i) + ".txt");
		vector<double> rhs2;
		shared_ptr<MatrixSkyline> matrix2{ new MatrixSkyline() };
		mesh3.refine_p();
		cout << "refine:\t" << i << endl;
		cout << "Number of the nodes:\t" << mesh3.GetNumberOfNodes() << endl;
		vector<double> solution2{ mesh3.getSolution() };
		FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline> method2{ problem.get(), &mesh3, matrix2.get(), &rhs2 };
		method2.Discretization();
		ESolver esl2{ Solvers::GMRES };
		solution2 = esl2.Solve(*matrix2.get(), rhs2, solution2, res, 100000, 1e-13);
		meths.push_back(method2);
		meshes.push_back(mesh3);
		FEAnalysis<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, CTriangularMesh, CTriangularMesh> anal;
		cout << method2.GetSolution(mesh3, solution2, Point(0.4, 0.4, 0.)) << endl;
		if (i > 0)
		{
			conv.push_back(anal.L2Norm(meths[4 + i], method2, meshes[4 + i], mesh3, all_souls[4 + i], solution2));
			cout << "Relative convergence in L2-norm:\t" << conv[4 + i] << endl;
			cout << "h-p convergence in L2-norm:\t" << anal.L2Norm(meths[4], method2, meshes[4], mesh3, all_souls[4], solution2) << endl;
			cout << "Convergence rate:\t" << (fabs(conv[3 + i] / conv[4 + i])) << "\t not log" << endl;
			cout << "Convergence rate:\t" << log2(fabs(conv[3 + i] / conv[4 + i])) << endl;
		}
		all_souls.push_back(solution2);
	}
	return 0;
}

const int test_case_elliptic_fem::elliptic_fem_solver() const
{
	cout << "Performing the test \"2d elliptic via fem solver\"" << endl;
	solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>> fem;
	CTriangularMesh mesh{ "grids//2d//unit_square_layer//unit.msh" };
	shared_ptr<CDiffusionScalar> problem{ new CDiffusionScalar() };
	problem->add_parameter(Terms::IDUDV, 16, 1);
	problem->add_parameter(Terms::IDUDV, 19, 1.5);
	problem->add_parameter(Terms::IDUDV, 20, 2);
	problem->add_parameter(Terms::IDUDV, 22, 1);
	problem->add_boundary_parameter(1, 1, 10);
	problem->add_boundary_parameter(1, 5, 20);
	vector<double> solution;
	fem.elliptic_solver(problem.get(), &mesh, &solution);
	cout << fem.get_value(mesh, solution, Point(0.5, 0.5, 0)) << endl;
	cout << fem.get_gradvalue(mesh, solution, Point(0.5, 0.5, 0)).x << endl;
	cout << fem.get_gradvalue(mesh, solution, Point(0.5, 0.5, 0)).y << endl;
	return 0;
}

const int test_case_elliptic_fem::elliptic_2layer_fem_2d_tria_h() const
{
	//solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>> fem;
	cout << "Performing the test \"2d elliptic via fem on triangles\"" << endl;
	CTriangularMesh mesh2{ "grids//2d//unit_incl//unit.msh" };
	shared_ptr<CDiffusionScalar> problem{ new CDiffusionScalar() };
	vector<double> res;
	vector<double> dg_solution;
	vector<vector<double>> all_souls;
	vector<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>> meths;
	vector<CTriangularMesh> meshes;
	vector<double> conv;
	const double T1 = 10;
	const double T2 = 20;
	const double T = 1. / ((T2 - T1) * (T2 - T1));
	problem->add_parameter(Terms::IDUDV, 18, 10.5);
	problem->add_parameter(Terms::IDUDV, 17, 0.5);
	problem->add_parameter(Terms::IDUDV, 16, 10.5);
	problem->add_boundary_parameter(1, 1, T1);
	problem->add_boundary_parameter(1, 3, T2);
	//if (problem->load_parameters("grids//2d//unit_square_layer//param.txt"))
		//return 1;
	{
		shared_ptr<MatrixSkyline> matrix{ new MatrixSkyline() };
		vector<double> rhs;
		CTriangularMesh mesh{ "grids//2d//unit_incl//unit.msh" };
		vector<double> solution{ mesh.getSolution() };
		FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline> method{ problem.get(), &mesh, matrix.get(), &rhs };
		method.Discretization();
		ESolver esl{ Solvers::GMRES };
		solution = esl.Solve(*matrix.get(), rhs, solution, res, 100000, 1e-13);
		all_souls.push_back(solution);
		meths.push_back(method);
		meshes.push_back(mesh);
		//fem.elliptic_solver(problem.get(), &mesh, &solution);
	}

	cout << "-----------------------------" << endl;
	cout << "h-refinement test convergence" << endl;
	cout << "-----------------------------" << endl;
	conv.push_back(0.);

	for (int i = 0; i < 4; ++i)
	{
		vector<double> rhs2;
		shared_ptr<MatrixSkyline> matrix2{ new MatrixSkyline() };
		mesh2.refine_h();
		cout << "refine:\t" << i << endl;
		cout << "Number of the nodes:\t" << mesh2.GetNumberOfNodes() << endl;
		vector<double> solution2{ mesh2.getSolution() };
		FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline> method2{ problem.get(), &mesh2, matrix2.get(), &rhs2 };
		method2.Discretization();
		ESolver esl2{ Solvers::GMRES };
		solution2 = esl2.Solve(*matrix2.get(), rhs2, solution2, res, 100000, 1e-13);
		meths.push_back(method2);
		meshes.push_back(mesh2);
		FEAnalysis<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, CTriangularMesh, CTriangularMesh> anal;

		conv.push_back(anal.L2Norm(meths[i], method2, meshes[i], mesh2, all_souls[i], solution2));
		cout << "Relative convergence in L2-norm:\t" << conv[i + 1] << endl;
		//cout << "Absolute convergence in L2-norm:\t" << anal.L2Norm(meths[0], method2, meshes[0], mesh2, all_souls[0], solution2) << endl;
		cout << "Convergence rate:\t" << (fabs(conv[i] / conv[i + 1])) << "\t not log" << endl;
		cout << "Convergence rate:\t" << log2(fabs(conv[i] / conv[i + 1])) << endl;
		all_souls.push_back(solution2);
	}
	double t2 = 0;
	cout << "-----------------------------" << endl;
	for (int i = 0; i < 4; ++i)
	{
		FEAnalysis<FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, FEMethod<CDiffusionScalar, CTriangularMesh, MatrixSkyline>, CTriangularMesh, CTriangularMesh> anal;
		const auto temp = anal.L2Norm(meths[i], meths[3], meshes[i], meshes[3], all_souls[i], all_souls[3]);
		cout << "Absolute convergence in L2-norm:\t" << temp << endl;
		cout << "Convergence rate:\t" << log2(fabs(t2 / temp)) << endl;
		t2 = temp;
	}
	return 0;
}

const int test_case_elliptic_fem::elliptic_fem_2d_rect_source() const
{
	auto bottom_f = [=](const int i, const Point& p)
	{
		//return 10;
		return p.x * p.x * p.x * p.x * p.x;
	};
	auto top_f = [=](const int i, const Point& p)
	{
		//return 20;
		return p.x * p.x * p.x * p.x * p.x;
	};
	auto fx = [=](const Point& p)
	{
		//return -2;
		return -20 * p.x * p.x * p.x;
	};
	cout << "Performing the test \"2d elliptic via fem solver\"" << endl;
	solvers::fem_solver<CDiffusionScalar, CRegularMesh, vector<double>> fem;
    const size_t nx = 3;
    const size_t ny = 3;
    CRegularMesh mesh{ 0, 0, 1, 1, int(std::pow(2, nx)), int(std::pow(2, nx)) };
	shared_ptr<CDiffusionScalar> problem{ new CDiffusionScalar() };
	auto source_f = [=](const int n, const Point& p)
	{
		double sum = 0.;
		const auto& elem = mesh.GetElement(n);
		const auto& dofs = elem->GetDoFs();
		for (size_t i = 0; i < dofs; ++i)
			sum += fx(mesh.GetNode(elem->GetNode(i))) * elem->GetShapeFunction(i, p);
		return sum;
	};
	parameter<double> source{ source_f };
	parameter<double> bottom{ bottom_f };
	parameter<double> top{ top_f };
	problem->addTerm(Terms::EFV);
	problem->add_parameter(Terms::IDUDV, 0, 1);
	problem->add_parameter(Terms::EFV, 0, source);
	problem->add_boundary_parameter(1, 0, bottom);
	problem->add_boundary_parameter(1, 1, top);
	problem->add_boundary_parameter(1, 2, top);
	problem->add_boundary_parameter(1, 3, top);
	vector<double> solution;
	fem.elliptic_solver(problem.get(), &mesh, &solution);
	cout << bottom_f(0, Point(0.5, 0.5, 0)) << endl;
	cout << fem.get_value(mesh, solution, Point(0.5, 0.5, 0)) << endl;
	cout << fem.get_gradvalue(mesh, solution, Point(0.5, 0.5, 0)).x << endl;
	cout << fem.get_gradvalue(mesh, solution, Point(0.5, 0.5, 0)).y << endl;
	return 0;
}


const int test_case_elliptic_fem::elliptic_gaussian_triangle() const
{
	return 0;
}

const int test_case_elliptic_fem::mass_matrix_3rd_order() const
{
    cout << "Triangle Mass Matrix Test" << endl;
    cout << "First order: ";
    int nodes[10];
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    nodes[6] = 6;
    nodes[7] = 7;
    nodes[8] = 8;
    nodes[9] = 9;
    vector<Point> points(3);
    points[0] = Point(4, 4);
    points[1] = Point(9, 4);
    points[2] = Point(4, 6);
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    points.push_back(Point(4, 5));
    CTriangle rect{ nodes, 10 };
    CTriangleLagrangeBasis basis{ &points[0], 3 };
    vector<vector<double>> MassMatrix = { { 2,1,1 },{ 1,2,1 },{ 1,1,2 } };
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            MassMatrix[i][j] *= basis.GetMeasure() / 24.;
    MassMatrix[0][0] = basis.GetMeasure() * (81. / 56 - 162. / 42 + 117. / 30 - 36. / 20 + 4. / 12) / 4;
    MassMatrix[0][0] = basis.GetMeasure() * (9. / 5040 - 3. / 1260 - 3. / 1260 + 1. / 360) * 81. / 4;
    double expected[3][3];
    auto func_mass = [&](const Point& p) {return basis.GetShapeFunction(5, p) * basis.GetShapeFunction(8, p); };
    expected[0][0] = rect.Integrate(func_mass, points);
    if (fabs(expected[0][0] - MassMatrix[0][0]) > 1e-10)
        {
            cout << basis.GetMeasure() << endl;
            cout << expected[0][0] << "\t" << MassMatrix[0][0] << endl;
            return 1;
        }
    /*for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            auto func_mass = [&](const Point& p) {return basis.GetShapeFunction(i, p) * basis.GetShapeFunction(j, p); };
            expected[i][j] = rect.Integrate(func_mass, points);
            if (fabs(expected[i][j] - MassMatrix[i][j]) > 1e-10)
                return 1;
        }*/
    cout << "Success" << endl;
    return 0;
    cout << "Second order: ";
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    CTriangle tr2{ nodes, 6 };
    CTriangleBasis bs2{ &points[0], 2 };
    auto func_mass2 = [&](const Point& p) {return basis.GetShapeFunction(3, p) * basis.GetShapeFunction(3, p); };
    const double val = tr2.Integrate(func_mass2, points);
    const double act = fabs(basis.GetMeasure()) * 4. / 720.;

    cout << "Third order: ";
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    points.push_back(Point(4, 5));
    nodes[6] = 6;
    nodes[7] = 7;
    nodes[8] = 8;
    nodes[9] = 9;
    CTriangle tr3{ nodes, 10 };
    CTriangleBasis bs3{ &points[0], 3 };
    auto func_mass3 = [&](const Point& p) {return bs3.GetShapeFunction(9, p) * bs3.GetShapeFunction(9, p); };
    const double val3 = tr3.Integrate(func_mass3, points);
    //const double act3 = fabs(bs3.GetMeasure()) * (2. * 48. / 40320. - 2 * 36. / 40320.);
    const double act3 = fabs(bs3.GetMeasure()) * (1./5040.);
    return 0;
}

const int test_case_elliptic_fem::strees_matrix_3rd_order() const
{
    cout << "Triangle Mass Matrix Test" << endl;
    cout << "First order: ";
    int nodes[10];
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    nodes[6] = 6;
    nodes[7] = 7;
    nodes[8] = 8;
    nodes[9] = 9;
    vector<Point> points(3);
    points[0] = Point(4, 4);
    points[1] = Point(9, 4);
    points[2] = Point(4, 6);
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    points.push_back(Point(4, 5));
    CTriangle rect{ nodes, 10 };
    CTriangleLagrangeBasis basis{ &points[0], 3 };
    vector<vector<double>> MassMatrix = { { 2,1,1 },{ 1,2,1 },{ 1,1,2 } };
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            MassMatrix[i][j] *= basis.GetMeasure() / 24.;
    MassMatrix[0][0] = basis.GetMeasure() * (81. / 56 - 162. / 42 + 117. / 30 - 36. / 20 + 4. / 12) / 4;
    double dx = 81. / 4 * (36 * basis.GetAlpha(0, 1) * basis.GetAlpha(1, 1) / 360 + 18 * basis.GetAlpha(1, 1) * basis.GetAlpha(2, 1) / 360 - 6 * basis.GetAlpha(1, 1) * basis.GetAlpha(2, 1) / 120 - 6 * basis.GetAlpha(0, 1) * basis.GetAlpha(1, 1) / 60 +
                           18 * basis.GetAlpha(0, 1) * basis.GetAlpha(2, 1) / 360 + 9 * basis.GetAlpha(2, 1) * basis.GetAlpha(2, 1) / 180 - 3 * basis.GetAlpha(2, 1) * basis.GetAlpha(2, 1) / 60 - 3 * basis.GetAlpha(0, 1) * basis.GetAlpha(2, 1) / 60 -
                           6 * basis.GetAlpha(0, 1) * basis.GetAlpha(1, 1) / 60 - 3 * basis.GetAlpha(1, 1) * basis.GetAlpha(2, 1) / 60 + basis.GetAlpha(1, 1) * basis.GetAlpha(2, 1) / 24 + basis.GetAlpha(0, 1) * basis.GetAlpha(1, 1) / 12 -
                           6 * basis.GetAlpha(0, 1) * basis.GetAlpha(2, 1) / 120 - 3 * basis.GetAlpha(2, 1) * basis.GetAlpha(2, 1) / 60 + basis.GetAlpha(2, 1) * basis.GetAlpha(2, 1) / 24 + basis.GetAlpha(0, 1) * basis.GetAlpha(2, 1) / 24);

    double dy = 81. / 4 * (36 * basis.GetAlpha(0, 2) * basis.GetAlpha(1, 2) / 360 + 18 * basis.GetAlpha(1, 2) * basis.GetAlpha(2, 2) / 360 - 6 * basis.GetAlpha(1, 2) * basis.GetAlpha(2, 2) / 120 - 6 * basis.GetAlpha(0, 2) * basis.GetAlpha(1, 2) / 60 +
                           18 * basis.GetAlpha(0, 2) * basis.GetAlpha(2, 2) / 360 + 9 * basis.GetAlpha(2, 2) * basis.GetAlpha(2, 2) / 180 - 3 * basis.GetAlpha(2, 2) * basis.GetAlpha(2, 2) / 60 - 3 * basis.GetAlpha(0, 2) * basis.GetAlpha(2, 2) / 60 -
                           6 * basis.GetAlpha(0, 2) * basis.GetAlpha(1, 2) / 60 - 3 * basis.GetAlpha(1, 2) * basis.GetAlpha(2, 2) / 60 + basis.GetAlpha(1, 2) * basis.GetAlpha(2, 2) / 24 + basis.GetAlpha(0, 2) * basis.GetAlpha(1, 2) / 12 -
                           6 * basis.GetAlpha(0, 2) * basis.GetAlpha(2, 2) / 120 - 3 * basis.GetAlpha(2, 2) * basis.GetAlpha(2, 2) / 60 + basis.GetAlpha(2, 2) * basis.GetAlpha(2, 2) / 24 + basis.GetAlpha(0, 2) * basis.GetAlpha(2, 2) / 24);
    double expected[3][3];
    MassMatrix[0][0] = basis.GetMeasure() * (dx + dy);
    auto func_stress = [&](const Point& p) {return basis.GetGradShapeFunction(5, p) * basis.GetGradShapeFunction(8, p); };
    expected[0][0] = rect.Integrate(func_stress, points);
    if (fabs(expected[0][0] - MassMatrix[0][0]) > 1e-10)
        {
            cout << basis.GetMeasure() << endl;
            cout << expected[0][0] << "\t" << MassMatrix[0][0] << endl;
            return 1;
        }
    /*for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            auto func_mass = [&](const Point& p) {return basis.GetShapeFunction(i, p) * basis.GetShapeFunction(j, p); };
            expected[i][j] = rect.Integrate(func_mass, points);
            if (fabs(expected[i][j] - MassMatrix[i][j]) > 1e-10)
                return 1;
        }*/
    cout << "Success" << endl;
    return 0;
    cout << "Second order: ";
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    nodes[3] = 3;
    nodes[4] = 4;
    nodes[5] = 5;
    CTriangle tr2{ nodes, 6 };
    CTriangleBasis bs2{ &points[0], 2 };
    auto func_mass2 = [&](const Point& p) {return basis.GetShapeFunction(3, p) * basis.GetShapeFunction(3, p); };
    const double val = tr2.Integrate(func_mass2, points);
    const double act = fabs(basis.GetMeasure()) * 4. / 720.;

    cout << "Third order: ";
    points.push_back(Point(6.5, 4));
    points.push_back(Point(6.5, 5));
    points.push_back(Point(4, 5));
    points.push_back(Point(4, 5));
    nodes[6] = 6;
    nodes[7] = 7;
    nodes[8] = 8;
    nodes[9] = 9;
    CTriangle tr3{ nodes, 10 };
    CTriangleBasis bs3{ &points[0], 3 };
    auto func_mass3 = [&](const Point& p) {return bs3.GetShapeFunction(9, p) * bs3.GetShapeFunction(9, p); };
    const double val3 = tr3.Integrate(func_mass3, points);
    //const double act3 = fabs(bs3.GetMeasure()) * (2. * 48. / 40320. - 2 * 36. / 40320.);
    const double act3 = fabs(bs3.GetMeasure()) * (1./5040.);
    return 0;
}

const int test_case_elliptic_fem::mass_matrix_4th_order() const
{
    return 0;
}

const int test_case_elliptic_fem::stress_matrix_4th_order() const
{
    return 1;
}


const int test_case_elliptic_fem::homotopy_conv_diff_fem(const double step) const
{
    ofstream ofs1("log.txt");
    std::random_device rd;
    std::mt19937 mt(rd());
    //std::uniform_real_distribution<double> distribution(, -4);
    const double a_max = 1;
    const int N = 10;
    std::normal_distribution<double> distribution(0, a_max);
    std::uniform_real_distribution<double> distribution_uni(0, 1);
    std::vector<double> sum(2);
    vector<double> temp(N);
    vector<Point> centers(N);
    ofs1 << "Generated Random Values" << endl;
    for (auto i = 0; i < N; ++i)
    {
        temp[i] = distribution(mt);
        centers[i].x = distribution_uni(mt);
        centers[i].y = distribution_uni(mt);
        ofs1 << temp[i] << endl;
    }
    ofs1 << "Centers X" << endl;
    for (auto i = 0; i < N; ++i)
        ofs1 << centers[i].x << endl;
    ofs1 << "Centers Y" << endl;
    for (auto i = 0; i < N; ++i)
        ofs1 << centers[i].y << endl;
    GaussianKernel gk{ N, centers };
    for (auto i = 0; i < 10; ++i)
        for (auto j = 0; j < 10; ++j)
        {
            cout << exp(gk.get_gp(temp, Point((double)i / 10., (double)j / 10.))) << endl;
        }
    cout.precision(15);
    scientific(cout);
    ofs1.precision(15);
    scientific(ofs1);
    //CTriangularMesh mesh{ "grids//2d//unit_square//mesh//0.msh" };
    CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2 };
    vector<vector<CTriangularMesh>> meshes;
    vector<vector<vector<double>>> solutions;
    vector<vector<double>> norms;
    vector<vector<double>> abs_norms;
    vector<vector<solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
    const int h_ref_max = 1;
    const int p_ref_max = 1;
    meshes.resize(p_ref_max);
    for (auto& it : meshes)
        it.resize(h_ref_max);
    fems.resize(p_ref_max);
    solutions.resize(p_ref_max);
    norms.resize(p_ref_max);
    abs_norms.resize(p_ref_max);
    for (auto i = 0; i < p_ref_max; ++i)
    {
        fems[i].resize(h_ref_max);
        solutions[i].resize(h_ref_max);
        norms[i].resize(h_ref_max);
        abs_norms[i].resize(h_ref_max);
        for (size_t j = 0; j < h_ref_max; ++j)
        {
            //meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
            meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j + 3), (int)std::pow(2, j + 3) };
            for (auto ii = 0; ii < i; ++ii)
                meshes[i][j].refine_p();
        }
    }
    double t = 0.;
    int nt = 1. / step;
    //while (fabs(t - (1. + step)) > 1e-13)
    for (int ii = 0; ii <= nt; ++ii)
    {
        ofs1 << "\nt:\t" << t << endl;
        for (auto i = 0; i < p_ref_max; ++i)
        {
            for (auto j = 0; j < h_ref_max; ++j)
            {
                abs_norms[i][j] = 0;
                norms[i][j] = 0;
            }
        }
        const double lal = 1e1;
        const double v = 1e2;
        const auto f0 = [=](const Point& p)
        {
            return -exp(p.x * p.y) * (p.x * p.x + p.y * p.y);
            return 0.;
        };
        const auto f1 = [=](const Point& p)
        {
            return v * (p.y * exp(p.x * p.y) + p.x * exp(p.x * p.y)) + lal * (p.y * p.y + p.x * p.x) * exp(p.x * p.y);
            //return 2. * v;
        };
        const auto fexp = [=](const Point& p)
        {
            return (1 - t) * f0(p) + t * f1(p);
            //return v * (p.y * exp(p.x * p.y) + p.x * exp(p.x * p.y)) + lal * (p.y * p.y + p.x * p.x) * exp(p.x * p.y);
            //return v;
            //return 20 * p.x + 2 * p.y + 0.4;
            //return 100.;
            //return 100 * exp(-10 * ((p.x - 0.5) * (p.x - 0.5) + (p.y - 0.5) * (p.y - 0.5)));
            //return std::pow((p.x + p.y - 10.), 4);
            return -3120. * std::pow((p.x + p.y - 10.), 38);
            //return -exp((p.x - 1) * (p.y - 1)) * ((p.x - 1) * (p.x - 1) + (p.y - 1) * (p.y - 1));
            //return -4.*(p.x - 10)*(p.y - 10) - (p.x - 10)*(2.*p.x + 2.*p.y - 20.) - (p.y - 10.)*(2.*p.x + 2.*p.y - 20.);
            //return -10.*std::pow(p.x+p.y-10.,4);
            //return -2.*exp(10. - p.x - p.y);
            //return 80.*exp(10 - p.x - p.y)*std::pow(p.x + p.y - 10., 39) - 3120.*exp(10. - p.y - p.x)*std::pow(p.x + p.y - 10., 38);
            //return exp(1 - p.x - p.y)*exp((p.x - 1) * (p.y - 1)) * (p.x - 1) + exp(1. - p.x - p.y) * exp((p.x - 1.) * (p.y - 1)) * (p.y - 1) - exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.x - 1) * (p.x - 1) -
            //	exp(1 - p.x - p.y) * exp((p.x - 1) * (p.y - 1)) * (p.y - 1) * (p.y - 1);
        };
        const auto bdn = [=](const int el, const int node, const Point& p)
        {
            return exp(p.x * p.y);
            return p.x + p.y;
            return p.x * p.x + p.y * p.y;
            //return exp((1 - p.x)*(1 - p.y));
            //return std::pow(10. - p.x - p.y, 40);
            //return 0.;
            //return std::pow(10. - p.x - p.y, 2);
            //return 10. - p.x - p.y;
        };
        const auto bdn1 = [=](const int el, const int node, const Point& p)
        {
            //return exp((1 - p.x)*(1 - p.y));
            return std::pow(10. - p.x - p.y, 40);
            //return 0.;
            //return std::pow(10. - p.x - p.y, 2);
            //return 10. - p.x - p.y;
        };
        const auto src = [=](const int el, const int node, const Point& p)
        {
            return fexp(p);
        };

        double _max = 0;
        double _min = 1000;
        const auto lam = [&](const int el, const int node, const Point& p)
        {
            //return exp(1. - p.x - p.y);
            //return std::pow(10. - p.x - p.y, 40);
            return lal;
            const double val = exp(gk.get_gp(temp, p));
            //const double val = gk.get_gp(temp, p);
            if (_max < val)
            {
                _max = val;
                cout << "max: " << val << endl;
            }
            if (_min > val)
            {
                _min = val;
                cout << "min: " << val << endl;
            }
            return val;
        };
        const auto vel = [&](const int el, const int node, const Point& p)
        {
            return t * Point(v, v);
        };
        const parameter<double> boundary_lin(bdn);
        const parameter<double> boundary_lin1(bdn1);
        const parameter<double> source(src);
        const parameter<double> lamda(lam);
        const parameter<Point>	velocity(vel);
        CDiffusionScalar problem;
        problem.addTerm(Terms::EFV);
        problem.add_parameter(Terms::EFV, 0, source);
        problem.addTerm(Terms::IDUV);
        problem.add_parameter(Terms::IDUV, 0, velocity);
        problem.add_parameter(Terms::IDUDV, 0, lamda);
        problem.add_boundary_parameter(1, 0, boundary_lin);
        problem.add_boundary_parameter(1, 1, boundary_lin);
        problem.add_boundary_parameter(1, 2, boundary_lin);
        problem.add_boundary_parameter(1, 3, boundary_lin);
        //fill(solutions[0][0].begin(), solutions[0][0].end(), 0);
        for (auto i = 0; i < p_ref_max; ++i)
        {
            for (auto j = 0; j < h_ref_max; ++j)
            {
                fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
                const auto& N_el = meshes[i][j].GetNumberOfElements();
                for (auto ii = 0; ii < N_el; ++ii)
                {
                    const auto& elem = meshes[i][j].GetElement(ii);
                    vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                    norms[i][j] += elem->Integrate([&](const Point& p)
                        {
                            const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p); return t * t;
                        }, pts);
                }
                norms[i][j] = sqrt(norms[i][j]);
                ofs1 << i << j << ": \t" << norms[i][j] << endl;
            }
        }
        const int Nx = 100;
        const int Ny = 100;
        for (auto i = 0; i < p_ref_max; ++i)
        {
            for (auto j = 0; j < h_ref_max; ++j)
            {
                ofstream res("results//mesh" + to_string(i) + to_string(j) + ".txt");
                for (auto ii = 0; ii < Nx + 1; ++ii)
                {
                    for (auto jj = 0; jj < Ny + 1; ++jj)
                    {
                        const auto xx = (double)ii / Nx;
                        const auto yy = (double)jj / Ny;
                        res << xx << "\t" << yy << "\t" << fems[i][j].get_value(meshes[i][j], solutions[i][j], Point((double)ii / Nx, (double)jj / Ny)) << endl;
                    }
                }
                //const auto& N_el = meshes[i][j].GetNumberOfElements();
            }
        }
        ofs1 << "Relative L2-norm" << endl;
        for (auto i = 0; i < p_ref_max; ++i)
        {
            for (auto j = 0; j < h_ref_max; ++j)
            {
                const auto& N_el = meshes[i][j].GetNumberOfElements();
                for (auto ii = 0; ii < N_el; ++ii)
                {
                    const auto& elem = meshes[i][j].GetElement(ii);
                    vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                    abs_norms[i][j] += elem->Integrate([&](const Point& p) {
                        const double t =
                            fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
                            //fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
                            bdn(ii, 0, p);
                        return t * t;
                        }, pts);
                }
                abs_norms[i][j] = sqrt(abs_norms[i][j]);
                ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
                //ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
            }
        }

        for (auto i = p_ref_max - 1; i < p_ref_max; ++i)
        {
            for (auto j = 0; j < h_ref_max - 1; ++j)
            {
                const auto& N_el = meshes[i][j].GetNumberOfElements();
                for (auto ii = 0; ii < N_el; ++ii)
                {
                    const auto& elem = meshes[i][j].GetElement(ii);
                    vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                    abs_norms[i][j] += elem->Integrate([&](const Point& p) {
                        const double t =
                            fems[i][j].get_value(meshes[i][j], solutions[i][j], p) -
                            //fems[p_ref_max - 1][h_ref_max - 1].get_value(meshes[p_ref_max - 1][h_ref_max - 1], solutions[p_ref_max - 1][h_ref_max - 1], p);
                            bdn(ii, 0, p);
                        return t * t;
                        }, pts);
                }
                abs_norms[i][j] = sqrt(abs_norms[i][j]);
                ofs1 << i << j << ": \t" << abs_norms[i][j] << endl;
                ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
            }
        }
        t += step;
        if (t > 1)
            t = 1;
    }
    ofs1 << "Convergence order for h" << endl;
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max - 1; ++j)
        {
            ofs1 << i << j << ": \t" << std::log2(abs_norms[i][j] / abs_norms[i][j + 1]) << endl;
        }
    }
    ofs1 << "Convergence order for p" << endl;
    for (auto i = 0; i < h_ref_max; ++i)
    {
        for (auto j = 0; j < p_ref_max - 1; ++j)
        {
            ofs1 << j << i << ": \t" << std::log2(abs_norms[j][i] / abs_norms[j + 1][i]) << endl;
        }
    }
    ofs1 << t << endl;
    return 0;
}

