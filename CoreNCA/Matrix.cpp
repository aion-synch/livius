#include "Matrix.h"

Algebra::Matrix::Matrix(const unsigned int & size, const std::vector<std::set<unsigned int>>& nonzero) :
	m_size{ size }
{
	m_elem.resize(size);
	for (unsigned int i = 0; i < size; ++i)
		m_elem[i].resize(size);
}

Algebra::Matrix::~Matrix()
{
	if (m_elem.size() > 0)
	{
		for (unsigned int i = 0; i < m_size; ++i)
			std::vector<double>().swap(m_elem[i]);
		std::vector<std::vector<double>>().swap(m_elem);
	}
}

void Algebra::Matrix::NullRow(const int row)
{
	for (unsigned int i = 0; i < m_size; ++i)
		m_elem[row][i] = 0;
	m_elem[row][row] = 1;
}

void Algebra::Matrix::NullMatrix()
{
	for(unsigned int i = 0; i < m_size; ++i)
		std::fill(m_elem[i].begin(), m_elem[i].end(), 0);
}

Algebra::Matrix & Algebra::Matrix::operator=(const Matrix & matrix)
{
	m_size = matrix.m_size;
	m_elem = matrix.m_elem;
	return *this;
}


Algebra::Matrix::Matrix(const Matrix & matrix)
{
	m_size = matrix.m_size;
	m_elem = matrix.m_elem;
}

void Algebra::Matrix::Create(const unsigned int & size, const std::vector<std::set<unsigned int>>& nonzero)
{
	m_elem.resize(size);
	for (unsigned int i = 0; i < size; ++i)
		m_elem[i].resize(size);
	m_size = size;
}

void Algebra::Matrix::Create(const unsigned int & size)
{
	m_elem.resize(size);
	for (unsigned int i = 0; i < size; ++i)
		m_elem[i].resize(size);
	m_size = size;
}
