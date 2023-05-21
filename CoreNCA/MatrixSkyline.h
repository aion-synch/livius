#ifndef CORENC_ALGEBRA_MATRIXSKYLINE_H_
#define CORENC_ALGEBRA_MATRIXSKYLINE_H_
#include <set>
#include <vector>
#include <memory>
#include "Matrix.h"
#include "MatrixDiag.h"
namespace Algebra
{
	class ESolver;
	enum class Solvers
	{
		BiCGStab,
		GMRES,
		GMRES_BiCGStab,
		Gauss,
		PARDISO
	};
    /**
     * The sparse (skyline) matrix format
     */
	class MatrixSkyline
	{
	public:
		MatrixSkyline(const unsigned int& Size, const std::vector<std::set<unsigned int>>& nonzero);
		MatrixSkyline(){};
		~MatrixSkyline();
		void				NullRow(int row);
		double& operator()(const int i, const int j)
		{
			int ind;
			/*for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
				if (m_colind[ind] == j)
					break;
			return (*this).m_valL[ind];*/
			if (i == j)
			{
				return (*this).m_valDiag[i];
			}
			if (i < j)
			{
				for (ind = m_rowptr[j]; ind < m_rowptr[j + 1]; ++ind)
					if (m_colind[ind] == i)
						break;
                return (*this).m_valU[ind];
			}
			else
			{
				for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
					if (m_colind[ind] == j)
						break;
				return (*this).m_valL[ind];
			}
			//return (*this)[];
		}
        const double operator()(const int i, const int j) const
        {
            int ind;
            /*for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
                if (m_colind[ind] == j)
                    break;
            return (*this).m_valL[ind];*/
            if (i == j)
            {
                return (*this).m_valDiag[i];
            }
            if (i < j)
            {
                for (ind = m_rowptr[j]; ind < m_rowptr[j + 1]; ++ind)
                    if (m_colind[ind] == i)
                        break;
                if (ind < m_rowptr[j + 1])
                    return (*this).m_valU[ind];
                return 0;
            }
            else
            {
                for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
                    if (m_colind[ind] == j)
                        break;
                if (ind < m_rowptr[i + 1])
                    return (*this).m_valL[ind];
                return 0;
            }
            return 0;
        }
        const double        GetElement(const int i, const int j) const
        {
            int ind;
            /*for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
                if (m_colind[ind] == j)
                    break;
            return (*this).m_valL[ind];*/
            if (i == j)
            {
                return (*this).m_valDiag[i];
            }
            if (i < j)
            {
                for (ind = m_rowptr[j]; ind < m_rowptr[j + 1]; ++ind)
                    if (m_colind[ind] == i)
                        break;
                if (ind < m_rowptr[j + 1])
                    return (*this).m_valU[ind];
                return 0;
            }
            else
            {
                for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
                    if (m_colind[ind] == j)
                        break;
                if (ind < m_rowptr[i + 1])
                    return (*this).m_valL[ind];
                return 0;
            }
            return 0;
        }
		const int			GetSize() const { return m_size; }
		void				NullMatrix();
		MatrixSkyline&		operator=(const MatrixSkyline&);
        //MatrixSkyline&		operator-(const MatrixSkyline&);
        //friend const MatrixSkyline operator-(const MatrixSkyline&, const MatrixSkyline&);
        // A - scal * B
        static const MatrixSkyline diff_skymatrix(const MatrixSkyline& matrix, const MatrixSkyline& B, const double scal)
        {
            MatrixSkyline C(matrix);
            if ((B.m_gsize == matrix.m_gsize) && (matrix.m_size == B.m_size))
            {
                for (int i = 0; i < B.m_gsize; ++i)
                {
                    C.m_valL[i] = matrix.m_valL[i] - scal * B.m_valL[i];
                    C.m_valU[i] = matrix.m_valU[i] - scal * B.m_valU[i];
                }
                for (int i = 0; i < B.m_size; ++i)
                {
                    C.m_valDiag[i] = matrix.m_valDiag[i] - scal * B.m_valDiag[i];
                }
            }
            return C;
        }
        static const MatrixSkyline diff_skymatrix(const MatrixSkyline& matrix, const MatrixSkyline& B, const double a, const double b)
        {
            MatrixSkyline C(matrix);
            if ((B.m_gsize == matrix.m_gsize) && (matrix.m_size == B.m_size))
            {
                for (int i = 0; i < B.m_gsize; ++i)
                {
                    C.m_valL[i] = a * matrix.m_valL[i] - b * B.m_valL[i];
                    C.m_valU[i] = a * matrix.m_valU[i] - b * B.m_valU[i];
                }
                for (int i = 0; i < B.m_size; ++i)
                {
                    C.m_valDiag[i] = a * matrix.m_valDiag[i] - b * B.m_valDiag[i];
                }
            }
            return C;
        }
        static const MatrixSkyline transpose_sky(const MatrixSkyline& matrix)
        {
            MatrixSkyline C(matrix);
            C.m_valL = matrix.m_valU;
            C.m_valU = matrix.m_valL;
            return C;
        }
		MatrixSkyline(const MatrixSkyline& matrix);
		void 				Create(const unsigned int& Size, const std::vector<std::set<unsigned int>>& nonzero);
		void				AddElement(const unsigned int i, const unsigned int j, const double a)
		{
			int ind;
			/*for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
			 if (m_colind[ind] == j)
			 break;
			 return (*this).m_valL[ind];*/
			if (i == j)
			{
				m_valDiag[i] += a;
				return;
			}
			if (i < j)
			{
				for (ind = m_rowptr[j]; ind < m_rowptr[j + 1]; ++ind)
					if (m_colind[ind] == i)
						break;
				m_valU[ind] += a;
				return;
			}
			else
			{
				for (ind = m_rowptr[i]; ind < m_rowptr[i + 1]; ++ind)
					if (m_colind[ind] == j)
						break;
				m_valL[ind] += a;
				return;
			}
		}
		
	private:
		//int*				m_rowptr = nullptr;
		std::vector<int>	m_rowptr;
		//int*				m_colind = nullptr;
		std::vector<int>	m_colind;
		//double*				m_valU = nullptr;
		std::vector<double>	m_valU;
		//double*				m_valL = nullptr;
		std::vector<double>	m_valL;
		//double*				m_valDiag = nullptr;
		std::vector<double>	m_valDiag;
		int					m_size;
		int					m_gsize;
		friend				ESolver;
		//friend				Matrix;
	};
	class ESolver
	{
	public:
		ESolver(const MatrixSkyline& matrix, const std::vector<double>& rightvector) :
			m_matrix{ matrix },
			m_rightvector( rightvector ),
			m_maxiter(20000),
			m_eps(1e-10)
		{
			m_solution.resize(matrix.m_size);
		}
		ESolver(){};
		ESolver(Solvers kek) :m_solver(kek){};
		void					Reload(const MatrixSkyline& matrix, const std::vector<double>& right);
		void					Solve(Solvers);
		const std::vector<double>		Solve(MatrixSkyline&, const std::vector<double>& rhs, std::vector<double>& sol, std::vector<double>& residual, const int iter, const double eps);
		const std::vector<double>		Solve(MatrixDiag&, const std::vector<double>& rhs, std::vector<double>& sol, std::vector<double>& residual, const int iter, const double eps);
		double					BiCGStab(const int _maxiter);
		double					GMRES(const int _maxiter);
		void					GMRES(MatrixSkyline&, const std::vector<double>& rhs, std::vector<double>& sol, std::vector<double>& residual, const int iter, const double eps);
		void					BiCGStab(MatrixSkyline&, const std::vector<double>& rhs, std::vector<double>& sol, std::vector<double>& residual, const int iter, const double eps);
		void					Gauss(MatrixSkyline&, const std::vector<double>& rhs, std::vector<double>& sol, std::vector<double>& residual, const int iter, const double eps);
		void					Gauss(Matrix&, const std::vector<double>& rhs, std::vector<double>& sol);
		void					Gauss(const Matrix&, double* in_out);
		void					Gauss(const Matrix&, double* in, double* out);
		void					Gauss(const Matrix&, const double* in, double* out);
		void					Pardiso(MatrixSkyline&, const std::vector<double>& rhs, std::vector<double>& sol);
		void					BiCGStabPrecond();
		const std::vector<double>	GetSolution() const{ return m_solution; }
		void					GetSolution(std::vector<double>& sol) const;
		void					MatrixprodVector(double*res, std::vector<double>& x, MatrixSkyline& m);
        void					MatrixprodVector(double*res, double* x, MatrixSkyline& m);
		void					MatrixprodVector(double*res, double* x, const Matrix& m);
		void					MatrixprodVector(double*res, const double* x, const Matrix& m);
		~ESolver();
	private:
		MatrixSkyline			m_matrix;
		std::vector<std::vector<double>> H, V, W;
		std::vector<double>		m_solution;
		std::vector<double>		m_rightvector;
		void					GMRESBiCGStab();
		void					MatrixprodVector(double*res, std::vector<double>& x, double*valDiag, double*valL, double*valU, int*rowptr, int*colind, int size);
		void					MatrixprodVector(double*res, std::vector<double>& x,
												std::vector<double>&valDiag, 
												std::vector<double>&valL,
												std::vector<double>&valU,
												std::vector<int>&rowptr,
												std::vector<int>&colind, int size);
		//void					MatrixprodVector(double*res, std::vector<double>& x, MatrixSkyline& m);
		void					MatrixprodVector(double*res, double* x,
												std::vector<double>&valDiag,
												std::vector<double>&valL,
												std::vector<double>&valU,
												std::vector<int>&rowptr,
												std::vector<int>&colind, int size);
		void					MatrixprodVector(double*res, double* x, double*valDiag, double*valL, double*valU, int*rowptr, int*colind, int size);
		void					MatrixprodVector(double* res, double* x, double* val, int* rowptr, int* colind, int size);
		const double			DotProd(double*, double*, int);
		void					zero_GMRES(std::vector<std::vector<double>>&, const int str, const int stl);
		void					mult_Ht_H_slae(double, double*, int);
		void					gauss(std::vector<std::vector<double>>&, double*, double*, int);
		int						find_max(std::vector<std::vector<double>>&, int, int);
		void					Copy(double *x, double *y, int n);
		void					mult_Vy(double*, double*, int);
		void					mult_Vy(double*, double*, int, int);
		const double			DotProd(const std::vector<double>&, const std::vector<double>&, int);
		const double			Norm(double*, int);
		const double			Norm(const std::vector<double>&, int);
		double					m_eps;
		Solvers					m_solver;
		int						m_maxiter;
		double*					m_LUvalL = nullptr;
		double*					m_LUvalU = nullptr;
		double*					m_LUvalDiag = nullptr;
		void					LUPrec();
		void					LSolve(double*, double*);
		void					USolve(double*, double*);
	public:
		auto					GetSolution() -> decltype(m_solution) { return m_solution; }
	};
}
#endif // CORENC_ALGEBRA_MATRIXSKYLINE_H_
