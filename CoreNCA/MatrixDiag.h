#ifndef CORENC_ALGEBRA_MATRIXDIAG_H
#define CORENC_ALGEBRA_MATRIXDIAG_H
#include <set>
#include <vector>

namespace Algebra
{
	class ESolver;
    /**
     * The diagonal matrix class
     */
	class MatrixDiag
	{
	public:
		MatrixDiag(const unsigned int& size, const std::vector<std::set<unsigned int>>& nonzero);
		MatrixDiag() {};
		~MatrixDiag();
		void									NullRow(const int row);
		double& operator()(const int i, const int j)
		{
			return (*this).m_valDiag[i];
		}
		const int								GetSize() const { return m_size; }
		void									NullMatrix();
		MatrixDiag&								operator=(const MatrixDiag&);
		MatrixDiag(const MatrixDiag& matrix);
		void 									Create(const unsigned int& size, const std::vector<std::set<unsigned int>>& nonzero);
		void									AddElement(const unsigned int i, const unsigned int j, const double a)
		{
			if (i == j)
			{
				m_valDiag[i] += a;
				return;
			}
			return;
		}
	private:
		std::vector<double>						m_valDiag;
		unsigned int							m_size{ 0 };
		friend									ESolver;
	};
}

#endif // !CORENC_ALGEBRA_MATRIXDIAG_H
