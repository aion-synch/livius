#include "test_case_solver.h"
#include "../CoreNCFEM/Grids/TriangularMesh.h"
#include "../CoreNCFEM/Grids/RegularMesh.h"
#include "../CoreNCFEM/Methods/FEMethod.h"
#include "../Problems/DiffusionScalar.h"
#include "../CoreNCA/MatrixSkyline.h"
#include "../CoreNCFEM/Methods/FEAnalysis.h"
#include "../Solvers/fem_solver.h"
#include "../CoreNCFEM/GaussianField.h"
//#include <eigen3/Eigen/Core>
//#include <Spectra/SymEigsSolver.h>
//#include <Spectra/SymGEigsSolver.h>
//#include <Spectra/GenEigsBase.h>

#include <random>
#define _USE_MATH_DEFINES

#include <math.h>
using namespace corenc;
using namespace std;
using namespace Mesh;
using namespace Algebra;
using namespace method;
//using namespace Spectra;

// M = diag(1, 2, ..., 10)
test_case_solver::test_case_solver()
{

}

test_case_solver::~test_case_solver()
{

}

/*const int avec_t(int n, double* a, double* b)
{
    // 2 x 2
    //  2 -1
    // -1  2
    b[0] = 2 * a[0] - a[1];
    b[1] = -a[0] + 2 * a[1];
    return 0;
}*/
const int solver(const Algebra::Matrix& matrix, double* x, double* res)
{
	Algebra::ESolver esl{ Algebra::Solvers::Gauss };
	esl.Gauss(matrix, x, res);
	return 0;
}

const int test_case_solver::gauss_solver() const
{
	cout << "Performing the test \"2d elliptic via fem Gauss solver\"" << endl;

      Matrix matrix;
      const int sz = 4;
      matrix.Create(sz);
      matrix(0, 0) = 1;
      matrix(0, 1) = 2;
      matrix(0, 2) = 3;
      matrix(0, 3) = 4;

      matrix(1, 0) = 5;
      matrix(1, 1) = 6;
      matrix(1, 2) = 7;
      matrix(1, 3) = 8;

      matrix(2, 0) = 9;
      matrix(2, 1) = 10;
      matrix(2, 2) = 11;
      matrix(2, 3) = 12;

      matrix(3, 0) = 13;
      matrix(3, 1) = 14;
      matrix(3, 2) = 15;
      matrix(3, 3) = 16;

      double x[sz];
      double res[sz];
      double tr[sz];
      x[0] = 1;
      x[1] = 1;
      x[2] = 1;
      x[3] = 1;
      solver(matrix, &x[0], &res[0]);
      
      for(int i = 0; i < sz; ++i)
      {
            for(int j = 0; j < sz; ++j)
                  tr[i] += matrix(i, j) * res[j];
            cout << res[i] << endl;
            if(fabs(tr[i] - x[i]) > 1e-13)
                  cout << "warning" << endl;
      }
	//dtest();
	//eigen_nsym();
	/*solvers::fem_solver<CDiffusionScalar, CTriangularMesh, vector<double>> fem;
	CTriangularMesh mesh{ Point{0,0}, Point{1,1}, 4, 4 };
	cout << mesh.GetNumberOfElements() << endl;
	shared_ptr<CDiffusionScalar> problem{ new CDiffusionScalar() };
	problem->add_parameter(Terms::IDUDV, 6, 1);
	problem->add_boundary_parameter(1, 1, 10);
	problem->add_boundary_parameter(1, 3, 20);
	vector<double> solution;
	fem.elliptic_solver(problem.get(), &mesh, &solution);
	cout << fem.get_value(mesh, solution, Point(0.5, 0.5, 0)) << endl;
	cout << fem.get_gradvalue(mesh, solution, Point(0.5, 0.5, 0)).x << endl;
	cout << fem.get_gradvalue(mesh, solution, Point(0.5, 0.5, 0)).y << endl;*/
	return 0;
}
