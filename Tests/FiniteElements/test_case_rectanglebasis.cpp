#include "test_case_rectanglebasis.h"
#include "../../CoreNCFEM/FiniteElements/Rectangle.h"
using namespace std;
using namespace corenc;
using namespace Mesh;
using namespace tests;


test_case_rectanglebasis::test_case_rectanglebasis()
{
}

test_case_rectanglebasis::~test_case_rectanglebasis()
{
}

const int test_case_rectanglebasis::mass_matrix() const
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
	vector<vector<double>> MassMatrix = { {4,2,2,1},{2,4,1,2},{2,1,4,2},{1,2,2,4} };
	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
			MassMatrix[i][j] *= hx * hy / 36.;
	CRectangleBasis basis{ &points[0], 4 };
	CRectangle rect{ nodes,4 };
	double expected[4][4];
	for(size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
		{
			auto func_mass = [&](const Point& p) {return basis.GetShapeFunction(i, p) * basis.GetShapeFunction(j, p); };
			expected[i][j] = rect.Integrate(func_mass, points);
			if (fabs(expected[i][j] - MassMatrix[i][j]) > 1e-13)
				return 1;
		}
	return 0;
}

const int test_case_rectanglebasis::stress_matrix() const
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
	CRectangleBasis basis{ &points[0], 4 };
	CRectangle rect{ nodes,4 };
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
