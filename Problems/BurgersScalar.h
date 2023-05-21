#ifndef CORENC_PROBLEMS_BURGERS_H_
#define CORENC_PROBLEMS_BURGERS_H_

#include "Problems.h"
#include <vector>
namespace corenc
{
	class CBurgersScalar : public CProblem
	{
	public:
		CBurgersScalar();
		~CBurgersScalar();
		Terms						getTerm(const unsigned int) const;
		const unsigned int			getNumberOfTerms() const;
		const int					setTerm(const unsigned int, const Terms&);
		const int					addTerm(const Terms&);
		const double				getFlux(const double) const;
		const int					removeTerm(const Terms&);
		const int					load_parameters(const std::string& file_name);
	private:
		std::vector<Terms>			m_terms;
	};
}
#endif // !CORENC_PROBLEMS_BURGERS_H_
