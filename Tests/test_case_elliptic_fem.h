#ifndef CORENC_TEST_CASE_ELLIPTIC_FEM_H_
#define CORENC_TEST_CASE_ELLIPTIC_FEM_H_

// SOME TEST PROBLEMS FOR ELLIPTIC CASE WITH FEM && DG\
// 0th, 1st, 2nd order definitely maybe more high-order
// LAGRANGE && HIERARHICAL BASIS FUNCTIONS
// LATER MAYBE EVEN TESTS WITH MULTISCALE

namespace corenc
{
	class test_case_elliptic_fem
	{
	public:
		test_case_elliptic_fem();
		~test_case_elliptic_fem();
		//const int					test_case_elliptic_fem_3d_tetra() const;
		const int					elliptic_fem_2d_tria() const;
		const int					elliptic_fem_solver() const;
		const int					elliptic_fem_square_lin_basis() const;
		const int					elliptic_fem_hp_fixed(const int h_ref_max, const int p_ref_max) const;
		const int					elliptic_fem_hp_fixed_triangle(const int h_ref_max, const int p_ref_max) const;
        const int					elliptic_fem_hp_lagrange_triangle(const int h_ref_max, const int p_ref_max) const;
		const int					elliptic_fem_hxhy_fixed_triangle(const int hx_max, const int hy_max) const;
		const int					conv_diff_fem_fixed_triangle(const int h_ref_max, const int p_ref_max) const;
        const int                   global_matrix(const int h_ref_max, const int p_ref_max) const;
		//const int					test_case_elliptic_fem_square_2nd_basis() const;
		//const int					test_case_elliptic_fem_square_nth_basis() const;
		const int					elliptic_2layer_fem_2d_tria_h() const;
		const int					elliptic_fem_2d_rect_source() const;
		const int					elliptic_gaussian_triangle() const;
        const int                   mass_matrix_3rd_order() const;
        const int                   strees_matrix_3rd_order() const;
        const int                   mass_matrix_4th_order() const;
        const int                   stress_matrix_4th_order() const;
        const int                   homotopy_conv_diff_fem(const double step) const;
		//const int					test_case_elliptic_fem_2d_rect() const;
		//const int					test_case_elliptic_fem_3d_hex() const;
		//const int					test_case_elliptic_dg_3d_tetra() const;
		//const int					test_case_elliptic_dg_2d_tria() const;
		//const int					test_case_elliptic_dg_2d_rect() const;
		//const int					test_case_elliptic_dg_3d_hex() const;
	};
}

#endif // !CORENC_TEST_CASE_ELLIPTIC_FEM_H_
