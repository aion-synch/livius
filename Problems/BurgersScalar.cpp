#include "BurgersScalar.h"
#include <vector>
using namespace corenc;

CBurgersScalar::CBurgersScalar()
{
	m_terms.resize(1);
	m_terms[0] = Terms::EUDV;
}

const double CBurgersScalar::getFlux(const double p) const
{
	return 0.5 * p * p;
}

const int CBurgersScalar::removeTerm(const Terms& term)
{
	for (unsigned i{ 0 }; i < m_terms.size(); ++i)
		if (m_terms[i] == term)
			m_terms.erase(m_terms.begin() + i);
	return 0;
}

Terms CBurgersScalar::getTerm(const unsigned int i) const
{
	if(i < m_terms.size())
		return m_terms[i];
	return m_terms[0];
}

const unsigned int CBurgersScalar::getNumberOfTerms() const
{
	return m_terms.size();
}

const int CBurgersScalar::setTerm(const unsigned int i, const Terms& term)
{
	if (i < m_terms.size())
	{
		m_terms[i] = term;
		return 0;
	}
	return 1;
}

const int CBurgersScalar::addTerm(const Terms& term)
{
	m_terms.push_back(term);
	return 0;
}

const int CBurgersScalar::load_parameters(const std::string& file_name)
{
	return 0;
}

CBurgersScalar::~CBurgersScalar()
{
	if (m_terms.size() > 0)
		std::vector<Terms>().swap(m_terms);
}
