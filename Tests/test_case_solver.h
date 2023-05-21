#ifndef CORENC_TEST_CASE_SOLVER_H_
#define CORENC_TEST_CASE_SOLVER_H_

// SOME TEST PROBLEMS FOR ELLIPTIC CASE WITH FEM && DG\
// 0th, 1st, 2nd order definitely maybe more high-order
// LAGRANGE && HIERARHICAL BASIS FUNCTIONS
// LATER MAYBE EVEN TESTS WITH MULTISCALE

namespace corenc
{
	class test_case_solver
	{
	public:
		test_case_solver();
		~test_case_solver();
		const int					gauss_solver() const;
	};
}

#endif // !CORENC_TEST_CASE_SOLVER_H_
