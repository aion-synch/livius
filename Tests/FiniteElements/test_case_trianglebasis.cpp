#include "test_case_trianglebasis.h"
#include "../../CoreNCFEM/FiniteElements/Triangle.h"
using namespace std;
using namespace corenc;
using namespace Mesh;
using namespace tests;


test_case_trianglebasis::test_case_trianglebasis()
{
}

test_case_trianglebasis::~test_case_trianglebasis()
{
}

const int test_case_trianglebasis::mass_matrix() const
{
	cout << "Triangle Mass Matrix Test" << endl;
	cout << "First order: ";
	int nodes[10];
	nodes[0] = 0;
	nodes[1] = 1;
	nodes[2] = 2;
	vector<Point> points(3);
	points[0] = Point(4, 4);
	points[1] = Point(9, 4);
	points[2] = Point(4, 6);
	CTriangle rect{ nodes,3 };
	CTriangleBasis basis{ &points[0], 1 };
	vector<vector<double>> MassMatrix = { { 2,1,1 },{ 1,2,1 },{ 1,1,2 } };
	for (size_t i = 0; i < 3; ++i)
		for (size_t j = 0; j < 3; ++j)
			MassMatrix[i][j] *= basis.GetMeasure() / 24.;
	double expected[3][3];
	for (size_t i = 0; i < 3; ++i)
		for (size_t j = 0; j < 3; ++j)
		{
			auto func_mass = [&](const Point& p) {return basis.GetShapeFunction(i, p) * basis.GetShapeFunction(j, p); };
			expected[i][j] = rect.Integrate(func_mass, points);
			if (fabs(expected[i][j] - MassMatrix[i][j]) > 1e-10)
				return 1;
		}
	cout << "Success" << endl;
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

const int test_case_trianglebasis::stress_matrix() const
{
	int nodes[4];
	nodes[0] = 0;
	nodes[1] = 1;
	nodes[2] = 2;
	nodes[3] = 3;
	vector<Point> points(4);
	points[0] = Point(4, 4);
	points[1] = Point(9, 4);
	points[2] = Point(4, 6);
	points[3] = Point(9, 6);
	const double hx = points[3].x - points[0].x;
	const double hy = points[3].y - points[0].y;
	vector<vector<double>> StressMatrix1 = { { 2,-2,1,-1 },{ -2,2,-1,1 },{ 1,-1,2,-2 },{ -1,1,-2,2 } };
	vector<vector<double>> StressMatrix2 = { { 2,1,-2,-1 },{ 1,2,-1,-2 },{ -2,-1,2,1 },{ -1,-2,1,2 } };
	vector<vector<double>> StressMatrix = { { 0,0,0,0 },{ 0,0,0,0 },{ 0,0,0,0 },{ 0,0,0,0 } };
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
			StressMatrix1[i][j] *= hy / hx / 6.;
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
			StressMatrix2[i][j] *= hx / hy / 6.;
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
			StressMatrix[i][j] = StressMatrix1[i][j] + StressMatrix2[i][j];
	CTriangleBasis basis{ &points[0], 4 };
	CTriangle rect{ nodes,4 };
	double expected[4][4];
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
		{
			auto func_stress = [&](const Point& p) {return basis.GetGradShapeFunction(i, p) * basis.GetGradShapeFunction(j, p); };
			expected[i][j] = rect.Integrate(func_stress, points);
			if (fabs(expected[i][j] - StressMatrix[i][j]) > 1e-13)
				return 1;
		}
	return 0;
}
