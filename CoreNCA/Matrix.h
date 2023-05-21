#ifndef CORENC_ALGEBRA_MATRIX_H
#define CORENC_ALGEBRA_MATRIX_H
#include <vector>
#include <set>
namespace Algebra
{
	class ESolver;
    /**
     * The Dense Matrix Class
     */
	class Matrix
	{
	public:
		Matrix(const unsigned int& size, const std::vector<std::set<unsigned int>>& nonzero);
		Matrix() {};
		~Matrix();
		void									NullRow(const int row);
		double& operator()(const int i, const int j)
		{
			return (*this).m_elem[i][j];
		}
		const int								GetSize() const { return m_size; }
		void									NullMatrix();
		Matrix&									operator=(const Matrix&);
		Matrix(const Matrix& matrix);
		void 									Create(const unsigned int& size, const std::vector<std::set<unsigned int>>& nonzero);
		void 									Create(const unsigned int& size);
		const double							GetElement(const int i, const int j)
		{
			return m_elem[i][j];
		}
		void									AddElement(const unsigned int i, const unsigned int j, const double a)
		{
			m_elem[i][j] += a;
			return;
		}
	private:
		std::vector<std::vector<double>>		m_elem;
		unsigned int							m_size{ 0 };
		friend									ESolver;
	};
}

#endif // !CORENC_ALGEBRA_MATRIX_H
