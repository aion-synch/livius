#pragma once
#ifndef CORENC_TEST_CASE_TRIANGLEBASIS_H_
#define CORENC_TEST_CASE_TRIANGLEBASIS_H_
namespace corenc
{
	namespace tests
	{
		class test_case_trianglebasis
		{
		public:
			test_case_trianglebasis();
			~test_case_trianglebasis();
			const int			mass_matrix() const;
			const int			stress_matrix() const;
		};
	}
}
#endif // !CORENC_TEST_CASE_TRIANGLEBASIS_H_
