#include "test_cases.h"
#include "test_case_elliptic_fem.h"
#include "test_case_solver.h"
#include "test_case_regular_mesh.h"
#include "FiniteElements/test_case_rectanglebasis.h"
#include "FiniteElements/test_case_trianglebasis.h"
#include <iostream>
#include <thread>
#include <future>
#include <chrono>
#include <ostream>
#include "../colors.h"
#include "test_conv_diff.h"
using namespace corenc;
using namespace std;
using namespace tests;

test_cases::test_cases()
{

}

test_cases::~test_cases()
{

}
const int test_cases::perform() const
{
    color::color_output("Performing tests...", color::YELLOW);
    test_case_elliptic_fem case1;
	test_case_regular_mesh case4;
	test_case_rectanglebasis case5;
    test_case_solver solver;
    test_conv_diff conv;
    test_case_solver gauss;
    //perform(std::bind(&test_case_elliptic_fem::elliptic_fem_2d_tria, std::ref(case1)));

    //case1.elliptic_fem_hp_fixed_triangle(4, 4);
    //case1.global_matrix(1, 4);
    //case1.strees_matrix_3rd_order();
    //case1.elliptic_fem_hp_lagrange_triangle(7, 1);
    //conv.conv_diff_fem(7);
    //case1.conv_diff_fem_fixed_triangle(7, 1);
    //case1.homotopy_conv_diff_fem(0.02);

    gauss.gauss_solver();
    perform(std::bind(&test_case_regular_mesh::construct_mesh, std::ref(case4)));
    perform(std::bind(&test_case_elliptic_fem::elliptic_fem_2d_rect_source, std::ref(case1)));
    perform(std::bind(&test_case_rectanglebasis::mass_matrix, std::ref(case5)));
	return 0;
}



const int test_cases::perform(const std::function<const int()>& f) const
{
	if (!f())
	{
		cout << "Success." << endl;
		return 0;
	}
	cout << "Something bad happened." << endl;
	return 1;
}

const int test_cases::perform(const std::function<const int(std::ostream& file_name)>& f, std::ostream& file_name) const
{
	if (!f(file_name))
	{
		cout << "Success." << endl;
		return 0;
	}
	cout << "Something bad happened." << endl;
	return 1;
}
