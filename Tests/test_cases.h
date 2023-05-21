#pragma once
#ifndef CORENC_TEST_CASES_H_
#define CORENC_TEST_CASES_H_
#include <functional>
#include <ostream>
namespace corenc
{
	class test_cases
	{
	public:
		test_cases();
		~test_cases();
        const int perform() const;
		const int perform(const std::function<const int()>&) const;
		const int perform(const std::function<const int(std::ostream&)>&, std::ostream&) const;
	};
}


#endif // !CORENC_TEST_CASES_H_
