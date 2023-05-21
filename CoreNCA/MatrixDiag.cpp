#include "MatrixDiag.h"

Algebra::MatrixDiag::MatrixDiag(const unsigned int & size, const std::vector<std::set<unsigned int>>& nonzero):
	m_size{size}
{
	m_valDiag.resize(size);
}

Algebra::MatrixDiag::~MatrixDiag()
{
	if (m_valDiag.size() > 0)
		std::vector<double>().swap(m_valDiag);
}

void Algebra::MatrixDiag::NullRow(const int row)
{
	m_valDiag[row] = 1;
}

void Algebra::MatrixDiag::NullMatrix()
{
	std::fill(m_valDiag.begin(), m_valDiag.end(), 0);
}

Algebra::MatrixDiag & Algebra::MatrixDiag::operator=(const MatrixDiag & matrix)
{
	m_size = matrix.m_size;
	m_valDiag = matrix.m_valDiag;
	return *this;
}


Algebra::MatrixDiag::MatrixDiag(const MatrixDiag & matrix)
{
	m_size = matrix.m_size;
	m_valDiag = matrix.m_valDiag;
}

void Algebra::MatrixDiag::Create(const unsigned int & size, const std::vector<std::set<unsigned int>>& nonzero)
{
	m_valDiag.resize(size);
	m_size = size;
}
