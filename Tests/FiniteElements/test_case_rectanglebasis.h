#pragma once
#ifndef CORENC_TEST_CASE_RECTANGLEBASIS_H_
#define CORENC_TEST_CASE_RECTANGLEBASIS_H_
namespace corenc
{
	namespace tests
	{
		class test_case_rectanglebasis
		{
		public:
			test_case_rectanglebasis();
			~test_case_rectanglebasis();
			const int			mass_matrix() const;
			const int			stress_matrix() const;
		};
	}
}
#endif // !CORENC_TEST_CASE_RECTANGLEBASIS_H_
