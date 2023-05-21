#pragma once
#ifndef CORENC_PROBLEMS_PROBLEMS_H_
#define CORENC_PROBLEMS_PROBLEMS_H_
#include "../CoreNCFEM/Point.h"
#include <string>

namespace corenc
{
	class CProblem
	{
	public:
		CProblem() {}
		virtual ~CProblem() {};
		virtual Terms					getTerm(const unsigned int) const = 0;
		virtual const unsigned int		getNumberOfTerms() const = 0;
		virtual const int				setTerm(const unsigned int, const Terms&) = 0;
		virtual const int				addTerm(const Terms&) = 0;
		virtual const int				load_parameters(const std::string& file_name) = 0;
	};
}

#endif // !CORENC_PROBLEMS_PROBLEMS_H_
