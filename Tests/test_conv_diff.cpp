#include "test_conv_diff.h"
#include "../CoreNCFEM/Grids/TriangularMesh.h"
#include "../CoreNCFEM/Grids/RegularMesh.h"
#include "../CoreNCFEM/Methods/FEMethod.h"
#include "../Problems/DiffusionScalar.h"
#include "../CoreNCA/MatrixSkyline.h"
#include "../CoreNCFEM/Methods/FEAnalysis.h"
#include "../Solvers/fem_solver.h"
#include "../Solvers/fem_solver_lib.h"
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

using namespace corenc;

void test_conv_diff::conv_diff_fem(const int h_ref_max, const int p_ref_max) const
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
    //CRegularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2 };
    CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 2, 2 };
    //vector<vector<CRegularMesh>> meshes;
    vector<vector<CTriangularMesh>> meshes;
    vector<vector<vector<double>>> solutions;
    vector<vector<double>> norms;
    vector<vector<double>> times;
    vector<vector<double>> abs_norms;
    //vector<vector<solvers::fem_solver<CDiffusionScalar, CRegularMesh, vector<double>>>> fems;
    vector<vector<solvers::fem_solver_lib<CDiffusionScalar, CTriangularMesh, vector<double>>>> fems;
    meshes.resize(p_ref_max);
    for (auto& it : meshes)
        it.resize(h_ref_max);
    fems.resize(p_ref_max);
    solutions.resize(p_ref_max);
    norms.resize(p_ref_max);
    times.resize(p_ref_max);
    abs_norms.resize(p_ref_max);
    for (auto i = 0; i < p_ref_max; ++i)
    {
        fems[i].resize(h_ref_max);
        solutions[i].resize(h_ref_max);
        norms[i].resize(h_ref_max);
        times[i].resize(h_ref_max);
        abs_norms[i].resize(h_ref_max);
        for (size_t j = 0; j < h_ref_max; ++j)
        {
            //meshes[i][j] = CTriangularMesh{ "grids//2d//unit_square//mesh//" + to_string(j) + ".msh" };
            //meshes[i][j] = CRegularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j + 1), (int)std::pow(2, j + 1) };
            meshes[i][j] = CTriangularMesh{ Point{0,0}, Point{1,1}, (int)std::pow(2, j + 1), (int)std::pow(2, j + 1) };
            for (auto ii = 0; ii < i; ++ii)
                meshes[i][j].refine_p();
        }
    }
    double lal = 1;
    lal = 1;
    double v = 2;
    double b = 2;
    v = 0;
    //v = 1024;
    //v = 1;
    const auto fexp = [=](const Point& p)
    {
        return 2 * b * b * lal * sin(b * M_PI * p.x) * M_PI * M_PI * sin(b * M_PI * p.y) + v * cos(b * M_PI * p.x) * b * M_PI * sin(b * M_PI * p.y) + v * cos(b * M_PI * p.y) * b * M_PI * sin(b * M_PI * p.x);
        //return 200 * lal * sin(10 * M_PI * p.x) * M_PI * M_PI * sin(10 * M_PI * p.y) + 10 * v * cos(M_PI * p.x) * M_PI * sin(M_PI * p.y) + 10 * v * sin(10 * M_PI * p.x) * cos(10 * M_PI * p.y);
        //return 0.;
        //return gk.get_gp(temp, Point(0.5, 0.5));
        auto um = exp(-(p.x - 0.5) * (p.x - 0.5) - (p.y - 0.5) * (p.y - 0.5));
        return v * (-2 * (p.x - 0.5) * um - 2 * (p.y - 0.5) * um) + lal * um * (-4 * p.x * p.x + 4 * p.x - 4 * p.y * p.y + 4 * p.y + 2);
        return -20 * v * um * ((p.x - 0.5) + (p.y - 0.5)) + lal * um * (-40 * p.x * p.x - 40 * p.x - 40 * p.y * p.y + 40 * p.y + 20);
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
        return sin(b * M_PI * p.x) * sin(b * M_PI * p.y);
        //return 10 * exp(-(p.x - 0.5) * (p.x - 0.5) - (p.y - 0.5) * (p.y - 0.5));
        //return 0.;
        //return 1.;
        return exp(-(p.x - 0.5) * (p.x - 0.5) - (p.y - 0.5) * (p.y - 0.5));
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
        //return 0.;
        return 2 * p.y;
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

    const auto gam = [=](const int el, const int node, const Point& p)
    {
        //diff
        //return -2.286578e1;
        //return -2.050554e1;
        //return -1.992979e1;
        //return -1.978679e1;
        // conv-diff
        //return -120.4649;
        //return -300.4829;
        //return -758.3057;
        return -677.4422;
        // h/4 h/8 h/16 h/32
        //return -7.949002e1;
        //return -1.271518e2;
        //return -1.201277e2;
        return -119.8063;
        //return -50.;
        //return -119.7392;
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
        //return Point(10, 10);
        return Point(v, v);
    };
    const parameter<double> boundary_lin(bdn);
    const parameter<double> boundary_lin1(1.);
    const parameter<double> boundary_lin2(1.);
    const parameter<double> boundary_lin3(bdn1);
    const parameter<double> source(src);
    const parameter<double> lamda(lam);
    const parameter<double> gamma(gam);
    const parameter<Point>	velocity(vel);
    CDiffusionScalar problem;
    problem.addTerm(Terms::EFV);
    //problem.addTerm(Terms::IUV);
    problem.add_parameter(Terms::EFV, 0, source);
    problem.addTerm(Terms::IDUV);
    //problem.addTerm(Terms::SUPG);
    problem.add_parameter(Terms::IDUV, 0, velocity);
    problem.add_parameter(Terms::IDUDV, 0, lamda);
    //problem.add_parameter(Terms::IUV, 0, gamma);
    problem.add_boundary_parameter(1, 0, boundary_lin);
    problem.add_boundary_parameter(1, 1, boundary_lin);
    problem.add_boundary_parameter(1, 2, boundary_lin);
    problem.add_boundary_parameter(1, 3, boundary_lin);

    for (auto i = 0; i < 1; ++i)
    {
        for (auto j = 0; j < 1; ++j)
        {
            chrono::steady_clock::time_point beg{ chrono::steady_clock::now() };
            fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
            chrono::steady_clock::time_point end{ chrono::steady_clock::now() };
            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg).count();
            times[i][j] = dur;
            const auto& N_el = meshes[i][j].GetNumberOfElements();
            for (auto ii = 0; ii < N_el; ++ii)
            {
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                // norm
                norms[i][j] += elem->Integrate([&](const Point& p)
                    {
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p, ii); return t * t;
                    }, pts);
                // diff norm
                /*norms[i][j] += elem->Integrate([&](const Point& p)
                    {
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p, ii) - bdn(0, 0, p); return t * t;
                    }, pts);*/
            }
            norms[i][j] = sqrt(norms[i][j]);
            ofs1 << i << j << ": \tPeclet:\t" << v * meshes[i][j].GetElement(0)->GetMeasure() / 2 / lal << "\tnomr\t" << norms[i][j] << endl;
            //ofs1 << i << j << ": \tPeclet:\t" << v * meshes[i][j].GetElement(0)->GetMeasure() / 6 << "\tnomr\t" << norms[i][j] << endl;
        }
    }
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 1; j < h_ref_max; ++j)
        {
            chrono::steady_clock::time_point beg{ chrono::steady_clock::now() };
            for (int k = 0; k < solutions[i][j].size(); ++k)
                solutions[i][j][k] = fems[i][j - 1].get_value(meshes[i][j - 1], solutions[i][j - 1], meshes[i][j].GetNode(k));
            fems[i][j].elliptic_solver(&problem, &meshes[i][j], &solutions[i][j]);
            chrono::steady_clock::time_point end{ chrono::steady_clock::now() };
            auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - beg).count();
            times[i][j] = dur;
            const auto& N_el = meshes[i][j].GetNumberOfElements();
            for (auto ii = 0; ii < N_el; ++ii)
            {
                const auto& elem = meshes[i][j].GetElement(ii);
                vector<Point> pts{ meshes[i][j].GetNode(elem->GetNode(0)), meshes[i][j].GetNode(elem->GetNode(1)), meshes[i][j].GetNode(elem->GetNode(2)) };
                // norm
                norms[i][j] += elem->Integrate([&](const Point& p)
                    {
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p, ii); return t * t;
                    }, pts);
                // diff norm
                /*norms[i][j] += elem->Integrate([&](const Point& p)
                    {
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p, ii) - bdn(0, 0, p); return t * t;
                    }, pts);*/
            }
            norms[i][j] = sqrt(norms[i][j]);
            ofs1 << i << j << ": \tPeclet:\t" << v * meshes[i][j].GetElement(0)->GetMeasure() / 2 / lal << "\tnomr\t" << norms[i][j] << endl;
            //ofs1 << i << j << ": \tPeclet:\t" << v * meshes[i][j].GetElement(0)->GetMeasure() / 6 << "\tnomr\t" << norms[i][j] << endl;
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
    ofs1 << endl;
    ofs1 << "times" << endl;
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max; ++j)
        {
            ofs1 << i << j << ": \t" << times[i][j] << endl;
            //ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
        }
    }
    ofs1 << endl;
    for (auto i = 0; i < p_ref_max; ++i)
    {
        for (auto j = 0; j < h_ref_max - 1; ++j)
        {
            ofs1 << i << j << ": \t" << log2(abs_norms[i][j] /  abs_norms[i][j + 1]) << endl;
            //ofs1 << i << j << ": \t" << abs_norms[i][j] / norms[p_ref_max - 1][h_ref_max - 1] << endl;
        }
    }
    const int Nx = 100;
    const int Ny = 100;
    /*for (auto i = 0; i < p_ref_max; ++i)
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
        }*/
    }


void test_conv_diff::conv_diff_eigen(const int h_ref_max, const int p_ref_max) const
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
        //return 1.198018e2;
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
    const auto gam = [=](const int el, const int node, const Point& p)
    {
        return 1.;
    };
    const parameter<double> boundary_lin(bdn);
    const parameter<double> boundary_lin1(bdn1);
    const parameter<double> source(src);
    const parameter<double> lamda(lam);
    const parameter<double> gamma(gam);
    const parameter<Point>	velocity(vel);
    CDiffusionScalar problem;
    problem.addTerm(Terms::EFV);
    problem.addTerm(Terms::IUV);
    problem.add_parameter(Terms::EFV, 0, source);
    problem.add_parameter(Terms::IUV, 0, gamma);
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
                        const double t = fems[i][j].get_value(meshes[i][j], solutions[i][j], p, ii); return t * t;
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
}
