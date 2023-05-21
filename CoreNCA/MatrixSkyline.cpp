#include "MatrixSkyline.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <ctime>
#include <cmath>
//#include <omp.h>
#define N_MIN 4096
#define _NOPE_
//#define _OMP_
//#define _MKL_
//#define _CUDA_
//#include <mkl.h>
#ifdef _CUDA_
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>
#endif // _CUDA_
using namespace Algebra;
using namespace std;
MatrixSkyline::MatrixSkyline(const unsigned int& N, const vector<set<unsigned int>>& a)
{
	//m_rowptr = new int[a.size() + 1];
	m_rowptr.resize(a.size() + 1);
	m_rowptr[0] = 0;
	/*for (int i = 0; i < a.size(); ++i)
	{
		m_rowptr[i + 1] = m_rowptr[i] + a[i].size();
	}
	m_colind = new int[m_rowptr[a.size()]];
	int k = 0;
	for (int i = 0; i < a.size(); ++i)
	{
		for (auto elem : a[i])
		{
			m_colind[k++] = elem;
		}
	}*/
	//m_rowptr[0] = 0;
	m_size = N;
	for (int i = 0; i < a.size(); ++i)
	{
		m_rowptr[i + 1] = m_rowptr[i] + a[i].size();
	}
	int size = m_rowptr[a.size()];
	//m_colind = new int[size];
	m_colind.resize(size);
	m_gsize = m_rowptr[a.size()];
	//m_valL = new double[m_gsize];
	m_valL.resize(m_gsize);
	for (int i = 0; i < m_gsize; ++i)
		m_valL[i] = 0;
	//cout << "size: " << m_size << endl;
	//cout << "nnz: " << m_gsize << endl;
	//m_valL = new double[size];
	//m_valU = new double[size];
	m_valU.resize(size);
	//m_valDiag = new double[N];
	m_valDiag.resize(N);
	for (unsigned int i = 0; i < N; ++i)
	{
		m_valDiag[i] = 0;
	}
	int k = 0;
	m_gsize = m_rowptr[a.size()];
	for (unsigned int i = 0; i < a.size(); ++i)
	{
		for (auto elem : a[i])
		{
			m_colind[k++] = elem;
		}
	}
	for (int i = 0; i < m_gsize; ++i)
		m_valU[i] = m_valL[i] = 0;
}

void MatrixSkyline::Create(const unsigned int& N, const vector<set<unsigned int>>& a)
{
	//m_rowptr = new int[a.size() + 1];
	m_rowptr.resize(a.size() + 1);
	m_rowptr[0] = 0;
	/*for (int i = 0; i < a.size(); ++i)
	 {
	 m_rowptr[i + 1] = m_rowptr[i] + a[i].size();
	 }
	 m_colind = new int[m_rowptr[a.size()]];
	 int k = 0;
	 for (int i = 0; i < a.size(); ++i)
	 {
	 for (auto elem : a[i])
	 {
	 m_colind[k++] = elem;
	 }
	 }*/
	//m_rowptr[0] = 0;
	m_size = N;
	for (int i = 0; i < a.size(); ++i)
	{
		m_rowptr[i + 1] = m_rowptr[i] + a[i].size();
	}
	int size = m_rowptr[a.size()];
	//m_colind = new int[size];
	m_colind.resize(size);
	m_gsize = m_rowptr[a.size()];
	//m_valL = new double[m_gsize];
	m_valL.resize(m_gsize);
	for (int i = 0; i < m_gsize; ++i)
		m_valL[i] = 0;
	//cout << "size: " << m_size << endl;
	//cout << "nnz: " << m_gsize << endl;
	//m_valL = new double[size];
	//m_valU = new double[size];
	m_valU.resize(size);
	//m_valDiag = new double[N];
	m_valDiag.resize(N);
	for (unsigned int i = 0; i < N; ++i)
	{
		m_valDiag[i] = 0;
	}
	int k = 0;
	m_gsize = m_rowptr[a.size()];
	for (unsigned int i = 0; i < a.size(); ++i)
	{
		for (auto elem : a[i])
		{
			m_colind[k++] = elem;
		}
	}
	for (int i = 0; i < m_gsize; ++i)
		m_valU[i] = m_valL[i] = 0;
}

MatrixSkyline::MatrixSkyline(const MatrixSkyline& matrix)
{
	m_gsize = matrix.m_gsize;
	m_size = matrix.m_size;
	m_rowptr = matrix.m_rowptr;
	m_colind = matrix.m_colind;
	m_valL = matrix.m_valL;
	m_valDiag = matrix.m_valDiag;
	m_valU = matrix.m_valU;
	/*m_rowptr = new int[m_size + 1];
	m_colind = new int[m_gsize];
	m_valL = new double[m_gsize];
	m_valU = new double[m_gsize];
	m_valDiag = new double[m_size];
	for (int i = 0; i < m_gsize; ++i)
	{
		m_valL[i] = matrix.m_valL[i];
		m_valU[i] = matrix.m_valU[i];
		m_colind[i] = matrix.m_colind[i];
	}
	for (int i = 0; i < m_size; ++i)
	{
		m_rowptr[i] = matrix.m_rowptr[i];
		m_valDiag[i] = matrix.m_valDiag[i];
	}
	m_rowptr[m_size] = matrix.m_rowptr[m_size];*/
}
void MatrixSkyline::NullRow(int row)
{
	int j = 0;
	/*for (int i = m_rowptr[row]; i < m_rowptr[row + 1]; ++i)
		if (m_colind[i] == row)
			m_valL[i] = 1;
		else
			m_valL[i] = 0;*/
	m_valDiag[row] = 1;
	for (; j < m_gsize; ++j)
	{
		if (m_colind[j] == row)
			m_valU[j] = 0;
	}
	for (j = m_rowptr[row]; j < m_rowptr[row + 1]; ++j)
		m_valL[j] = 0;
}
void MatrixSkyline::NullMatrix()
{
	for (int i = 0; i < m_gsize; ++i)
	{
		m_valL[i] = 0.;
		m_valU[i] = 0.;
	}
	for (int i = 0; i < m_size; ++i)
	{
		m_valDiag[i] = 0.;
	}
}
MatrixSkyline& MatrixSkyline::operator=(const MatrixSkyline& matrix)
{
    m_gsize = matrix.m_gsize;
    m_size = matrix.m_size;
    m_valL = matrix.m_valL;
    m_valU = matrix.m_valU;
    m_colind = matrix.m_colind;
    m_rowptr = matrix.m_rowptr;
    m_valDiag = matrix.m_valDiag;
    /*if ((m_gsize == matrix.m_gsize) && (matrix.m_size == m_size))
	{
		//m_valL = matrix.m_valL;
		//m_valU = matrix.m_valU;
		//m_colind = matrix.m_colind;
		//m_rowptr = matrix.m_rowptr;
		//m_valDiag = matrix.m_valDiag;
		for (int i = 0; i < m_gsize; ++i)
		{
			m_valL[i] = matrix.m_valL[i];
			m_valU[i] = matrix.m_valU[i];
			m_colind[i] = matrix.m_colind[i];
		}
		for (int i = 0; i < m_size; ++i)
		{
			m_rowptr[i] = matrix.m_rowptr[i];
			m_valDiag[i] = matrix.m_valDiag[i];
		}
		m_rowptr[m_size] = matrix.m_rowptr[m_size];
    }*/
	return *this;
}


MatrixSkyline::~MatrixSkyline()
{
	if (m_valU.size() > 0)
		std::vector<double>().swap(m_valU);
	if (m_valL.size() > 0)
		std::vector<double>().swap(m_valL);
	if (m_rowptr.size() > 0)
		std::vector<int>().swap(m_rowptr);
	if (m_valDiag.size() > 0)
		std::vector<double>().swap(m_valDiag);
	if (m_colind.size() > 0)
		std::vector<int>().swap(m_colind);
	/*if (m_valU)
		delete[] m_valU;
	if (m_valL)
		delete[] m_valL;
	if (m_rowptr)
		delete[] m_rowptr;
	if (m_colind)
		delete[] m_colind;
	if (m_valDiag)
		delete[] m_valDiag;
	m_valU = nullptr;
	m_valL = nullptr;
	m_rowptr = nullptr;
	m_colind = nullptr;
	m_valDiag = nullptr;*/
}
void ESolver::Reload(const MatrixSkyline& matrix, const vector<double>& right)
{
	//for (unsigned i = 0; i < right.size(); ++i)
		//m_rightvector[i] = right[i];
	m_rightvector = right;
	m_matrix = matrix;
}
void ESolver::GetSolution(vector<double>& vec) const
{
	vec.resize(m_solution.size());
	for (unsigned int i = 0; i < m_solution.size(); ++i)
		vec[i] = m_solution[i];
}
void ESolver::Solve(Solvers solver)
{
	switch (solver)
	{
	case Solvers::BiCGStab:
	{
		BiCGStab(m_maxiter);
		break;
	}
	case Solvers::GMRES:
	{
		GMRES(m_maxiter);
		break;
	}
	case Solvers::GMRES_BiCGStab:
	{
		GMRESBiCGStab();
		break;
	}
	default:
		break;
	}
}
void ESolver::GMRESBiCGStab()
{
	int iter{ 0 };
	int m_size = m_matrix.GetSize();
	double res;
	ofstream outfile{ "residual.txt" };
	while (iter < m_maxiter)
	{
		if ((res = fabs(GMRES(10))) < m_eps)
			break;
		outfile << iter << "\t" << res << endl;
		if ((res = fabs(BiCGStab(40))) < m_eps)
			break;
		outfile << iter << "\t" << res << endl;
		++iter;
	}
	outfile.close();
	cout << "Number of iteration:\t" << iter << endl;
}

double ESolver::GMRES(const int _maxiter)
{
	double *r, *y, *vec, *Res0, *Res, *x, betta, residual{1e10}, res_r, res_r_;
	int i, j, l, flag{ 0 }, flag_m, old_h_m, iter{ 0 }, h_m{30};
	old_h_m = 1;
	int m_size = m_matrix.GetSize();
	H.resize(h_m + 1);
	for (i = 0; i < h_m + 1; ++i)
		H[i].resize(h_m);
	V.resize(h_m);
	for (i = 0; i < h_m; ++i)
		V[i].resize(m_size);
	W.resize(h_m);
	for (i = 0; i < h_m; ++i)
		W[i].resize(m_size);
	r = new double[m_size];
	vec = new double[m_size];
	Res0 = new double[m_size];
	Res = new double[m_size];
	y = new double[m_size];
	x = new double[m_size];
	for (i = 0; i < m_size; ++i)
		Res0[i] = 0;
	res_r = residual / 10;
	res_r_ = 0;
	for (i = 0; i < m_size; ++i)
		Res[i] = m_solution[i];
	MatrixprodVector(r, Res, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
	for (i = 0; i < m_size; ++i)
		r[i] = m_rightvector[i] - r[i];
	//ofstream outfile{ "residual.txt" };
	while (flag == 0 && iter < _maxiter)
	{
		zero_GMRES(H, h_m + 1, h_m);
		zero_GMRES(V, h_m, m_size);
		zero_GMRES(W, h_m, m_size);

		flag_m = 0;
		h_m = old_h_m;
		betta = Norm(r, m_size);
		for (i = 0; i < m_size; ++i)
			V[0][i] = r[i] / betta;
		for (j = 0; j < h_m && flag_m == 0; ++j)
		{
			MatrixprodVector(&W[j][0], &V[j][0], m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
			for (i = 0; i <= j; ++i)
			{
				H[i][j] = DotProd(&W[j][0], &V[i][0], m_size);
				for (l = 0; l < m_size; ++l)
					W[j][l] = W[j][l] - H[i][j] * V[i][l];
			}
			H[j + 1][j] = Norm(&W[j][0], m_size);
			if (fabs(H[j + 1][j]) < 1e-10)
			{
				h_m = j;
				flag_m = 1;
			}
			else
				for (l = 0; l < m_size; ++l)
					V[j + 1][l] = W[j + 1][l] / H[j + 1][j];
		}
		mult_Ht_H_slae(betta, y, h_m);
		mult_Vy(y, vec, h_m);
		for (l = 0; l < m_size; ++l)
			Res[l] = Res[l] + vec[l];
		MatrixprodVector(r, Res, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
		for (i = 0; i < m_size; ++i)
			r[i] = m_rightvector[i] - r[i];
		residual = res_r;
		res_r = Norm(r, m_size);
		if (res_r < m_eps)
			flag = 1;
		else
			for (i = 0; i < m_size; ++i)
				m_solution[i] = Res[i];
		//outfile << iter << "\t" << res_r << endl;
		++iter;
	}
	//outfile.close();
	//for (i = 0; i < m_size; ++i)
	//	m_solution[i] = Res0[i];
	delete[] Res;
	delete[] r;
	delete[] vec;
	delete[] Res0;
	delete[] y;
	delete[] x;
	return res_r;
}
void ESolver::mult_Ht_H_slae(double betta, double *y, int m)
{
	int i, j, k;
	vector<vector<double>> HH;
	HH.resize(m);
	for (i = 0; i<m; i++) HH[i].resize(m);

	double *g, s;

	g = new double[m];

	for (i = 0; i<m; i++)
	{
		for (j = 0; j<m; j++)
		{
			s = 0;
			for (k = 0; k<m + 1; k++)	s += H[k][i] * H[k][j];
			HH[i][j] = s;
		}
		g[i] = H[0][i] * betta;
	}
	gauss(HH, y, g, m);
	delete[] g;
}
void ESolver::gauss(vector<vector<double>>&C, double *y, double *f, int m)
{
	int i, j, k;
	double *tmp, t, kf;
	tmp = new double[m];
	for (i = 0; i<m - 1; i++)
	{
		j = find_max(C, i, m);
		if (j != i)
		{
			Copy(tmp, &C[j][0], m);
			Copy(&C[j][0], &C[i][0], m);
			Copy(&C[i][0], tmp, m);
			t = f[i];
			f[i] = f[j];
			f[j] = t;
		}
		for (j = i + 1; j<m; j++)
		{
			kf = C[j][i] / C[i][i];
			C[j][i] = 0;
			for (k = i + 1; k<m; k++)
			{
				C[j][k] -= kf*C[i][k];
			}
			f[j] -= kf*f[i];
		}
	}
	for (i = m - 1; i >= 0; i--)
	{
		y[i] = f[i] / C[i][i];
		for (j = i - 1; j >= 0; j--)
			f[j] -= C[j][i] * y[i];
	}
	delete[] tmp;
}
int ESolver::find_max(vector<vector<double>>&C, int j, int m)
{
	int max, i;
	max = j;
	for (i = j + 1; i<m; i++)
		if (fabs(C[i][j])>fabs(C[max][j]))
			max = i;
	return max;
}
void ESolver::mult_Vy(double *y, double *g, int m)
{
	int i, j;
	double s;
	for (i = 0; i<m_matrix.m_size; i++)
	{
		s = 0;
		for (j = 0; j<m; j++)
			s += V[j][i] * y[j];
		g[i] = s;
	}
}
void ESolver::mult_Vy(double *y, double *g, int m, int size)
{
	int i, j;
	double s;
	for (i = 0; i<size; i++)
	{
		s = 0;
		for (j = 0; j<m; j++)
			s += V[j][i] * y[j];
		g[i] = s;
	}
}
void ESolver::zero_GMRES(vector<vector<double>>& B, const int str, const int stl)
{
	int i, j;
	for (i = 0; i<str; i++)
		for (j = 0; j<stl; j++)
			B[i][j] = 0;
}
void ESolver::Copy(double *x, double *y, int n)
{
	for (int i = 0; i<n; i++) x[i] = y[i];
}
double ESolver::BiCGStab(const int _maxiter)
{
	double *r, *rt, *s, *p, *v, *t;
	double alp, bet, ro, omega, omegal, rol;
	double norm;
	int k = 0;
	int m_size = m_matrix.GetSize();
	r = new double[m_size];
	rt = new double[m_size];
	p = new double[m_size];
	v = new double[m_size];
	s = new double[m_size];
	t = new double[m_size];
	int incx = 1;
	int incy = 1;
	double residual;
	char transa = 'N';
	double alpha = 1;
	double beta = 0;
	//char matdescraL = 'TLNF';
	//char matdescraU = 'TUUF';
	ro = alp = omega = omegal = rol = 1;
	// v = p = 0;
	double *temp = new double[m_size];
	m_solution.resize(m_size);
	int i = 0;
	//MatrixprodVector(temp, m_solution, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
#ifdef _CUDA_
	cudaError_t cudaStat1, cudaStat2, cudaStat3, cudaStat4, cudaStat5, cudaStat6; 
	cusparseStatus_t status; 
	cusparseHandle_t handle = 0; 
	cublasHandle_t bhandle = 0;
	cusparseMatDescr_t descr = 0;
	cublasCreate(&bhandle);
	cusparseCreate(&handle);
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
	double* valL;
	int* rowptr;
	double *temp2;
	int* colind;
	double* temp3;
	double* x;
	cudaStat1 = cudaMalloc(&valL, m_matrix.m_gsize*sizeof(double));
	cudaStat1 = cudaMalloc(&rowptr, (m_size+1)*sizeof(int));
	cudaStat1 = cudaMalloc(&colind, m_matrix.m_gsize*sizeof(int));
	cudaStat1 = cudaMalloc(&temp2, m_size*sizeof(double));
	cudaStat1 = cudaMalloc(&x, m_size*sizeof(double));
	cudaStat1 = cudaMalloc(&temp3, m_size*sizeof(double));
	cudaStat1 = cudaMemcpy(valL, m_matrix.m_valL, m_matrix.m_gsize*sizeof(double), cudaMemcpyHostToDevice);
	cudaStat1 = cudaMemcpy(rowptr, m_matrix.m_rowptr, (m_size+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaStat1 = cudaMemcpy(colind, m_matrix.m_colind, m_matrix.m_gsize*sizeof(int), cudaMemcpyHostToDevice);
	for (i = 0; i < m_size; ++i)
		temp[i] = 0;
	cudaStat1 = cudaMemcpy(temp2, temp, m_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaStat1 = cudaMemcpy(temp3, temp, m_size*sizeof(double), cudaMemcpyHostToDevice);
	int m_gsize = m_matrix.m_gsize;
	status = cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_size, 
		m_size, m_gsize, &alpha, descr, valL, rowptr, colind, temp2, &beta, temp3);
	cudaStat1 = cudaMemcpy(temp, temp3, m_size*sizeof(double), cudaMemcpyDeviceToHost);
	double* d_rt;
	double* d_r;
	double* d_v;
	double* d_p;
	double* d_s;
	double* d_t;
	double* d_rhs;
	double negomega;
	double negalp;
	double one = 1;
	double none = -1;
	cudaStat1 = cudaMalloc(&d_t, m_size*sizeof(double));
	cudaStat1 = cudaMalloc(&d_rt, m_size * sizeof(double));
	cudaStat1 = cudaMalloc(&d_r, m_size * sizeof(double));
	cudaStat1 = cudaMalloc(&d_p, m_size*sizeof(double));
	cudaStat1 = cudaMalloc(&d_v, m_size*sizeof(double));
	cudaStat1 = cudaMalloc(&d_s, m_size*sizeof(double));
	cudaStat1 = cudaMalloc(&d_rhs, m_size*sizeof(double));
	//cudaMemcpy(d_rt, rt, m_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rhs, &*m_rightvector.begin(), m_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_r, &*m_rightvector.begin(), m_size*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_r, r, m_size*sizeof(double), cudaMemcpyHostToDevice);
	cublasDaxpy(bhandle, m_size, &none, temp3, 1, d_r, 1);
	cudaMemcpy(d_rt, d_r, m_size*sizeof(double), cudaMemcpyDeviceToDevice);
	cublasDnrm2(bhandle, m_size, d_r, 1, &norm);
	cublasDnrm2(bhandle, m_size, d_rhs, 1, &residual);
	residual = norm / residual;
	//cudaMemcpy(d_v, v, m_size*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_p, p, m_size*sizeof(double), cudaMemcpyHostToDevice);
#endif // _CUDA_
#if defined(_NOPE_) || defined(_OMP_)
	//MatrixprodVector(temp, &*m_solution.begin(), m_matrix.m_valL, m_matrix.m_rowptr, m_matrix.m_colind, m_matrix.m_size);
	MatrixprodVector(temp, m_solution, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
#endif // _NOPE_
	//for (int i = 0; i < m_matrix.m_gsize; ++i)
		//++m_matrix.m_rowptr[i];
	//mkl_dskymv(&transa, &m_size, &m_size, &alpha, &matdescraL, m_matrix.m_valL, m_matrix.m_rowptr, &*m_solution.begin(), &beta, temp);
	//mkl_dskymv(&transa, &m_size, &m_size, &alpha, &matdescraU, m_matrix.m_valU, m_matrix.m_rowptr, &*m_solution.begin(), &beta, temp2);
	//for (int i = 0; i < m_size; ++i)
		//temp[i] += temp2[i];
	//for (int i = 0; i < m_matrix.m_gsize; ++i)
		//--m_matrix.m_rowptr[i];
#ifdef _MKL_
	mkl_cspblas_dcsrgemv(&transa, &m_size, m_matrix.m_valL, m_matrix.m_rowptr, m_matrix.m_colind, &*m_solution.begin(), temp);
#endif // _MKL_
#ifndef _CUDA_
	for (int i = 0; i < m_size; ++i)
	{
		rt[i] = r[i] = m_rightvector[i] - temp[i];
		p[i] = v[i] = 0;
	}
	residual = Norm(r, m_size) / Norm(m_rightvector, m_size);
	
	//for /(int i = 0; i < m_size; ++i)
		//cout << m_rightvector[i] << endl;
#endif // _CUDA_
	k = 1;
	//ofstream outfile{ "residual.txt" };
	
	unsigned int start_time = clock();
	//int m_endtime;
	while ((k < _maxiter) && (residual > m_eps))
	{
		rol = ro; omegal = omega;
#if defined(_NOPE_) || defined(_OMP_)
		ro = DotProd(rt, r, m_size);
#endif // _NOPE_
#ifdef _MKL_
		ro = ddot(&m_size, rt, &incx, r, &incx);
#endif // _MKL_
#ifdef _CUDA_
		cublasDdot(bhandle, m_size, d_rt, 1, d_r, 1, &ro);
#endif // _CUDA_
		if (ro == 0)
			cout << "FAIL!" << endl;
		if (k == 1)
		{
#ifndef _CUDA_
			for (i = 0; i < m_size; ++i)
				p[i] = r[i];
#endif // _CUDA_
#ifdef _CUDA_
			cublasDcopy(bhandle, m_size, d_r, 1, d_p, 1);
#endif // _CUDA_
		}
		else
		{
			bet = (ro / rol)*(alp / omegal);
#ifndef _CUDA_
			for (i = 0; i < m_size; ++i)
				p[i] = r[i] + bet * (p[i] - omegal * v[i]);
#endif // _CUDA_
#ifdef _CUDA_
			negomega = -omegal;
			cublasDaxpy(bhandle, m_size, &negomega, d_v, 1, d_p, 1);
			cublasDscal(bhandle, m_size, &bet, d_p, 1);
			cublasDaxpy(bhandle, m_size, &one, d_r, 1, d_p, 1);
#endif // _CUDA_
		}
		//mkl_dskymv(&transa, &m_size, &m_size, &alpha, &matdescraL, m_matrix.m_valL, m_matrix.m_colind, p, &beta, v);
		//MatrixprodVector(v, p, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
#ifdef _MKL_
		mkl_cspblas_dcsrgemv(&transa, &m_size, m_matrix.m_valL, m_matrix.m_rowptr, m_matrix.m_colind, p, v);
		alp = ro / ddot(&m_size, rt, &incx, v, &incy);
#endif // _MKL_
#if defined(_NOPE_) || defined(_OMP_)
		//MatrixprodVector(temp, m_solution, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
		MatrixprodVector(v, p, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
		//MatrixprodVector(v, p, m_matrix.m_valL, m_matrix.m_rowptr, m_matrix.m_colind, m_matrix.m_size);
		alp = ro / DotProd(rt, v, m_size);
#endif // _NOPE_
#ifdef _CUDA_
		cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_size,
			m_size, m_gsize, &alpha, descr, valL, rowptr, colind, d_p, &beta, d_v);	
		cublasDdot(bhandle, m_size, d_rt, 1, d_v, 1, &alp);
		alp = ro / alp;
		negalp = -alp;
		cublasDaxpy(bhandle, m_size, &negalp, d_v, 1, d_r, 1);
		cublasDaxpy(bhandle, m_size, &alp, d_p, 1, x, 1);
		cublasDnrm2(bhandle, m_size, d_r, 1, &norm);
#endif // _CUDA_
#ifndef _CUDA_
		for (i = 0; i < m_size; ++i)
			s[i] = r[i] - alp * v[i];
		norm = Norm(s, m_size);
#endif // _CUDA_
		if (norm < m_eps)
		{
#ifndef _CUDA_
			for (i = 0; i < m_size; ++i)
				m_solution[i] += alp*p[i];
#endif // _CUDA_
			break;
		}
		//MatrixprodVector(t, s, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
#ifdef _MKL_
		mkl_cspblas_dcsrgemv(&transa, &m_size, m_matrix.m_valL, m_matrix.m_rowptr, m_matrix.m_colind, s, t);
#endif // _MKL_
#if defined(_NOPE_) || defined(_OMP_)
		MatrixprodVector(t, s, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_size);
		//MatrixprodVector(t, s, m_matrix.m_valL, m_matrix.m_rowptr, m_matrix.m_colind, m_matrix.m_size);
		omega = DotProd(t, s, m_size) / DotProd(t, t, m_size);
#endif // _NOPE_
#ifdef _MKL_
		omega = ddot(&m_size, t, &incx, s, &incy) / ddot(&m_size, t, &incx, t, &incy);
#endif // _MKL_
#ifdef _CUDA_
		cusparseDcsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_size,
			m_size, m_gsize, &alpha, descr, valL, rowptr, colind, d_r, &beta, d_t);
		cublasDdot(bhandle, m_size, d_t, 1, d_r, 1, &temp[0]);
		cublasDdot(bhandle, m_size, d_t, 1, d_t, 1, &temp[1]);
		omega = temp[0] / temp[1];
		negomega = -omega;
		cublasDaxpy(bhandle, m_size, &omega, d_r, 1, x, 1);
		cublasDaxpy(bhandle, m_size, &negomega, d_t, 1, d_r, 1);
		cublasDnrm2(bhandle, m_size, d_r, 1, &norm);
		cublasDnrm2(bhandle, m_size, d_rhs, 1, &residual);
		residual = norm / residual;
#endif // _CUDA_
#ifndef _CUDA_
		for (i = 0; i < m_size; ++i)
		{
			m_solution[i] += alp * p[i] + omega * s[i];
			r[i] = s[i] - omega * t[i];
		}
		residual = Norm(r, m_size) / Norm(m_rightvector, m_size);
#endif // _CUDA_
		++k;
	}
#ifdef _CUDA_
	cudaMemcpy(&*m_solution.begin(), x, m_size*sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(x);
	cudaFree(d_v);
	cudaFree(d_p);
	cudaFree(d_s);
	cudaFree(d_t);
	cudaFree(valL);
	cudaFree(rowptr);
	cudaFree(colind);
	cudaFree(temp2);
	cudaFree(temp3);
	cusparseDestroyMatDescr(descr);
	descr = 0;
	cusparseDestroy(handle);
	cublasDestroy(bhandle);
	handle = 0;
#endif
	//cout << "residual: " << residual << endl;
	//outfile.close();
	delete[] temp;
	delete[] r;
	delete[] rt;
	delete[] s;
	delete[] t;
	delete[] v;
	delete[] p;
	return residual;
}
const double ESolver::Norm(double*u, int size)
{
	const int incx = 1;
#ifdef _MKL_
	return sqrt(ddot(&size, u, &incx, u, &incx));
#endif // _MKL_
	return sqrt(DotProd(u, u, size));
}
const double ESolver::Norm(const vector<double>& u, int size)
{
	const int incx = 1;
#ifdef _MKL_
	return sqrt(ddot(&size, &*u.begin(), &incx, &*u.begin(), &incx));
#endif // _MKL_
	return sqrt(DotProd(u, u, size));
}
const double ESolver::DotProd(double *u, double *v, int size)
{
	double res = 0;
	int i = 0;
#ifdef _NOPE_
	for (i = 0; i < size; ++i)
	{
		res += u[i] * v[i];
	}
#endif // _NOPE_
#ifdef _OMP_
	if (size > N_MIN)
	{
		#pragma omp parallel for shared(u, v) private(i) reduction(+:res)
		for (i = 0; i < size; ++i)
			res += u[i] * v[i];
	}
	else
	{
		for (i = 0; i < size; ++i)
			res += u[i] * v[i];
	}
#endif // _OMP_
	return res;
}
const double ESolver::DotProd(const vector<double>& u, const vector<double>& v, int size)
{
	double res = 0;
	int i = 0;
#ifdef _NOPE_
	for (i = 0; i < size; ++i)
	{
		res += u[i] * v[i];
	}
#endif // _NOPE_
#ifdef _OMP_
	if (size > N_MIN)
	{
#pragma omp parallel for shared(u, v) private(i) reduction(+:res)
		for (i = 0; i < size; ++i)
			res += u[i] * v[i];
	}
	else
	{
		for (i = 0; i < size; ++i)
			res += u[i] * v[i];
	}
#endif // _OMP_
	return res;
}
void ESolver::MatrixprodVector(double*res, vector<double>& x, double*valDiag, double*valL, double*valU, int*rowptr, int*colind, int size)
{
	int fst, lst, j;
	for (int i = 0; i < size; i++)
	{
		res[i] = valDiag[i] * x[i];
		fst = rowptr[i];
		lst = rowptr[i + 1];
		for (int k = fst; k < lst; k++)
		{
			j = colind[k];
			res[i] += valL[k] * x[j];
			res[j] += valU[k] * x[i];
		}
	}
}
void ESolver::MatrixprodVector(double*res, double* x,
	vector<double>&valDiag,
	vector<double>&valL,
	vector<double>&valU,
	vector<int>&rowptr,
	vector<int>&colind, int size)
{
	int fst, lst, j;
	for (int i = 0; i < size; i++)
	{
		res[i] = valDiag[i] * x[i];
		fst = rowptr[i];
		lst = rowptr[i + 1];
		for (int k = fst; k < lst; k++)
		{
			j = colind[k];
			res[i] += valL[k] * x[j];
			res[j] += valU[k] * x[i];
		}
	}
}
void ESolver::MatrixprodVector(double*res, vector<double>& x, MatrixSkyline& matrix)
{
	MatrixprodVector(res, x, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, matrix.m_size);
}

void ESolver::MatrixprodVector(double*res, double* x, MatrixSkyline& matrix)
{
    MatrixprodVector(res, x, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, matrix.m_size);
}

void ESolver::MatrixprodVector(double* res, double* x, const Matrix& matrix)
{
	//res = new double[matrix.m_size];
	for(int i = 0; i < matrix.m_size; ++i)
	{
		res[i] = 0;
		for(int j = 0; j < matrix.m_size; ++j)
			res[i] += matrix.m_elem[i][j] * x[j];
	}
}

void ESolver::MatrixprodVector(double* res, const double* x, const Matrix& matrix)
{
	//res = new double[matrix.m_size];
	for(int i = 0; i < matrix.m_size; ++i)
	{
		res[i] = 0;
		for(int j = 0; j < matrix.m_size; ++j)
			res[i] += matrix.m_elem[i][j] * x[j];
	}
}

void ESolver::MatrixprodVector(double*res, vector<double>& x,
	vector<double>&valDiag,
	vector<double>&valL,
	vector<double>&valU,
	vector<int>&rowptr,
	vector<int>&colind, int size)
{
	int fst, lst, j;
	for (int i = 0; i < size; i++)
	{
		res[i] = valDiag[i] * x[i];
		fst = rowptr[i];
		lst = rowptr[i + 1];
		for (int k = fst; k < lst; k++)
		{
			j = colind[k];
			res[i] += valL[k] * x[j];
			res[j] += valU[k] * x[i];
		}
	}
}

void ESolver::MatrixprodVector(double*res, double* x, double*valDiag, double*valL, double*valU, int*rowptr, int*colind, int size)
{
	int fst, lst, j;
	for (int i = 0; i < size; i++)
	{
		res[i] = valDiag[i] * x[i];
		fst = rowptr[i];
		lst = rowptr[i + 1];
		for (int k = fst; k < lst; k++)
		{
			j = colind[k];
			res[i] += valL[k] * x[j];
			res[j] += valU[k] * x[i];
		}
	}
}
void ESolver::MatrixprodVector(double* res, double* x, double* val, int* rowptr, int* colind, int size)
{
#ifndef _OMP_
	for (int i = 0; i < size; ++i)
	{
		res[i] = 0;
		for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
			res[i] += x[colind[j]] * val[j];
	}
#else
	int i, j;
#pragma omp parallel shared(rowptr, colind, val, x, res) private(i, j)
	{
#pragma omp for
		for (i = 0; i < size; ++i)
		{
			res[i] = 0;
			for (j = rowptr[i]; j < rowptr[i + 1]; ++j)
				res[i] += x[colind[j]] * val[j];
		}
	}

#endif // _OMP_

}
const std::vector<double> ESolver::Solve(MatrixSkyline& matrix, const std::vector<double>& rightvector, std::vector<double>& solution, std::vector<double>& res, const int maxiter, const double eps)
{
    if (solution.size() != matrix.GetSize())
        solution.resize(matrix.GetSize());
	m_solution.resize(matrix.GetSize());
    m_solution = solution;
	res.resize(matrix.GetSize());
	switch (m_solver)
	{
	case Solvers::BiCGStab:
	{
		BiCGStab(matrix, rightvector, solution, res, maxiter, eps);
		break;
	}
	case Solvers::GMRES:
	{
		GMRES(matrix, rightvector, solution, res, maxiter, eps);
		break;
	}
	case Solvers::Gauss:
	{
		Gauss(matrix, rightvector, solution, res, maxiter, eps);
		break;
	}
	case Solvers::PARDISO:
	{
		Pardiso(matrix, rightvector, solution);
	}
	}
	/*double *r, *rt, *s, *p, *v, *t;
	double alp, bet, ro, omega, omegal, rol;
	double norm;
	int k = 0;
	int m_size = matrix.GetSize();
	r = new double[m_size];
	rt = new double[m_size];
	p = new double[m_size];
	v = new double[m_size];
	s = new double[m_size];
	t = new double[m_size];
	int incx = 1;
	int incy = 1;
	double residual;
	char transa = 'N';
	double alpha = 1;
	double beta = 0;
	char matdescraL = 'TLNF';
	char matdescraU = 'TUUF';
	ro = alp = omega = omegal = rol = 1;
	// v = p = 0;
	double *temp = new double[m_size];
	m_solution.resize(m_size);
	m_solution = solution;
	int i = 0;
	MatrixprodVector(temp, m_solution, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
	for (int i = 0; i < m_size; ++i)
	{
		rt[i] = r[i] = rightvector[i] - temp[i];
		p[i] = v[i] = 0;
	}
	residual = Norm(r, m_size) / Norm(rightvector, m_size);

	k = 1;

	unsigned int start_time = clock();
	int m_endtime;
	while ((k < maxiter) && (residual > eps))
	{
		rol = ro; omegal = omega;
		ro = DotProd(rt, r, m_size);
		if (ro == 0)
			cout << "FAIL!" << endl;
		if (k == 1)
		{
			for (i = 0; i < m_size; ++i)
				p[i] = r[i];
		}
		else
		{
			bet = (ro / rol)*(alp / omegal);
			for (i = 0; i < m_size; ++i)
				p[i] = r[i] + bet * (p[i] - omegal * v[i]);
		}
		MatrixprodVector(v, p, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
		alp = ro / DotProd(rt, v, m_size);
		for (i = 0; i < m_size; ++i)
			s[i] = r[i] - alp * v[i];
		norm = Norm(s, m_size);
		if (norm < m_eps)
		{
			for (i = 0; i < m_size; ++i)
				m_solution[i] += alp*p[i];
			break;
		}
		MatrixprodVector(t, s, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
		omega = DotProd(t, s, m_size) / DotProd(t, t, m_size);
		for (i = 0; i < m_size; ++i)
		{
			m_solution[i] += alp * p[i] + omega * s[i];
			r[i] = s[i] - omega * t[i];
		}
		residual = Norm(r, m_size) / Norm(rightvector, m_size);
		++k;
	}
	cout << "residual: " << residual << endl;
	delete[] temp;
	delete[] r;
	delete[] rt;
	delete[] s;
	delete[] t;
	delete[] v;
	delete[] p;*/

	return m_solution;
}

void ESolver::GMRES(MatrixSkyline& matrix, const std::vector<double>& rightvector, std::vector<double>& solution, std::vector<double>& res, const int maxiter, const double eps)
{
	solution.resize(matrix.GetSize());
	m_solution.resize(matrix.GetSize());
	//m_solution = solution;
	res.resize(matrix.GetSize());
	double *r, *y, *vec, *Res0, *Res, *x, betta, residual{ 1e10 }, res_r, res_r_;
	int i, j, l, flag{ 0 }, flag_m, old_h_m, iter{ 0 }, h_m{ 30 };
	old_h_m = 1;
	int m_size = matrix.GetSize();
	H.resize(h_m + 1);
	for (i = 0; i < h_m + 1; ++i)
		H[i].resize(h_m);
	V.resize(h_m);
	for (i = 0; i < h_m; ++i)
		V[i].resize(m_size);
	W.resize(h_m);
	for (i = 0; i < h_m; ++i)
		W[i].resize(m_size);
	r = new double[m_size];
	vec = new double[m_size];
	Res0 = new double[m_size];
	Res = new double[m_size];
	y = new double[m_size];
	x = new double[m_size];
	for (i = 0; i < m_size; ++i)
		Res0[i] = 0;
	res_r = residual / 10;
	res_r_ = 0;
	for (i = 0; i < m_size; ++i)
		Res[i] = m_solution[i];
	MatrixprodVector(r, Res, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
	for (i = 0; i < m_size; ++i)
		r[i] = rightvector[i] - r[i];
	//ofstream outfile{ "residual.txt" };
	while (flag == 0 && iter < maxiter)
	{
		zero_GMRES(H, h_m + 1, h_m);
		zero_GMRES(V, h_m, m_size);
		zero_GMRES(W, h_m, m_size);

		flag_m = 0;
		h_m = old_h_m;
		betta = Norm(r, m_size);
		for (i = 0; i < m_size; ++i)
			V[0][i] = r[i] / betta;
		for (j = 0; j < h_m && flag_m == 0; ++j)
		{
			MatrixprodVector(&W[j][0], &V[j][0], matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
			for (i = 0; i <= j; ++i)
			{
				H[i][j] = DotProd(&W[j][0], &V[i][0], m_size);
				for (l = 0; l < m_size; ++l)
					W[j][l] = W[j][l] - H[i][j] * V[i][l];
			}
			H[j + 1][j] = Norm(&W[j][0], m_size);
			if (fabs(H[j + 1][j]) < 1e-10)
			{
				h_m = j;
				flag_m = 1;
			}
			else
				for (l = 0; l < m_size; ++l)
					V[j + 1][l] = W[j + 1][l] / H[j + 1][j];
		}
		mult_Ht_H_slae(betta, y, h_m);
		mult_Vy(y, vec, h_m, m_size);
		for (l = 0; l < m_size; ++l)
			Res[l] = Res[l] + vec[l];
		MatrixprodVector(r, Res, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
		for (i = 0; i < m_size; ++i)
			r[i] = rightvector[i] - r[i];
		residual = res_r;
		res_r = Norm(r, m_size);
		if (res_r < eps)
			flag = 1;
		else
			for (i = 0; i < m_size; ++i)
				m_solution[i] = Res[i];
		//outfile << iter << "\t" << res_r << endl;
		++iter;
	}
	if (iter >= maxiter)
	{
		cout << "WARNING!" << endl;
		cout << res_r << endl;
	}
	//outfile.close();
	//for (i = 0; i < m_size; ++i)
	//	m_solution[i] = Res0[i];
	MatrixprodVector(r, m_solution, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
	for (int i = 0; i < m_size; ++i)
		res[i] = rightvector[i] - r[i];
	//cout << "The residual: " << residual << "\t Iterations: " << iter << endl;
	delete[] Res;
	delete[] r;
	delete[] vec;
	delete[] Res0;
	delete[] y;
	delete[] x;
}

void ESolver::BiCGStab(MatrixSkyline& matrix, const std::vector<double>& rightvector, std::vector<double>& solution, std::vector<double>& res, const int maxiter, const double eps)
{
	{
		double *r, *rt, *s, *p, *v, *t;
		double alp, bet, ro, omega, omegal, rol;
		double norm;
		int k = 0;
		int m_size = matrix.GetSize();
		r = new double[m_size];
		rt = new double[m_size];
		p = new double[m_size];
		v = new double[m_size];
		s = new double[m_size];
		t = new double[m_size];
		int incx = 1;
		int incy = 1;
		double residual;
		char transa = 'N';
		double alpha = 1;
		double beta = 0;
		//char matdescraL = 'TLNF';
		//char matdescraU = 'TUUF';
		ro = alp = omega = omegal = rol = 1;
		// v = p = 0;
		double *temp = new double[m_size];
		int i = 0;
		MatrixprodVector(temp, m_solution, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
		for (int i = 0; i < m_size; ++i)
		{
			rt[i] = r[i] = rightvector[i] - temp[i];
			p[i] = v[i] = 0;
		}
		residual = Norm(r, m_size) / Norm(rightvector, m_size);

		k = 1;

		unsigned int start_time = clock();
		//int m_endtime;
		while ((k < maxiter) && (residual > eps))
		{
			rol = ro; omegal = omega;
			ro = DotProd(rt, r, m_size);
			if (ro == 0)
				cout << "FAIL!" << endl;
			if (k == 1)
			{
				for (i = 0; i < m_size; ++i)
					p[i] = r[i];
			}
			else
			{
				bet = (ro / rol)*(alp / omegal);
				for (i = 0; i < m_size; ++i)
					p[i] = r[i] + bet * (p[i] - omegal * v[i]);
			}
			MatrixprodVector(v, p, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
			alp = ro / DotProd(rt, v, m_size);
			for (i = 0; i < m_size; ++i)
				s[i] = r[i] - alp * v[i];
			norm = Norm(s, m_size);
			if (norm < m_eps)
			{
				for (i = 0; i < m_size; ++i)
					m_solution[i] += alp * p[i];
				break;
			}
			MatrixprodVector(t, s, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
			omega = DotProd(t, s, m_size) / DotProd(t, t, m_size);
			for (i = 0; i < m_size; ++i)
			{
				m_solution[i] += alp * p[i] + omega * s[i];
				r[i] = s[i] - omega * t[i];
			}
			residual = Norm(r, m_size) / Norm(rightvector, m_size);
            //cout << residual << endl;
			++k;
		}
		MatrixprodVector(temp, m_solution, matrix.m_valDiag, matrix.m_valL, matrix.m_valU, matrix.m_rowptr, matrix.m_colind, m_size);
		if (k >= maxiter)
		{
			cout << "WARNING!" << endl;
			cout << residual << endl;
		}
        //cout << "iter: " << k << endl;
		for (int i = 0; i < m_size; ++i)
			res[i] = rightvector[i] - temp[i];
        //cout << "The residual: " << residual <<"\t Iterations: " << k << endl;
		delete[] temp;
		delete[] r;
		delete[] rt;
		delete[] s;
		delete[] t;
		delete[] v;
		delete[] p;
	}
}
void ESolver::Gauss(MatrixSkyline& matrix, const std::vector<double>& rightvector, std::vector<double>& solution, std::vector<double>& res, const int maxiter, const double eps)
{
	cout << "Gauss solver is currently not implemented for skyline matrices!" << endl;
}

void ESolver::Pardiso(MatrixSkyline& matrix, const std::vector<double>& rhs, std::vector<double>& sol)
{
    //_MKL_DSS_HANDLE_t pt[64];
    //for (size_t i = 0; i < 64; ++i)
        //pt[i] = 0;
	int maxfct = 1;
	int mnum = 1;
	int mtype = 11;
	int phase = 13;
	vector<double> a;
	int ind = 0;

	//a.resize(m_valL.size() * 2 + m_valDiag.size());
	vector<int> ia;
	vector<int> ja;
	ia.push_back(1);
	for (int i = 0; i < matrix.m_size; ++i)
	{
		for (int j = 0; j < matrix.m_size; ++j)
		{
			if (i == j)
			{
				a.push_back(matrix.m_valDiag[i]);
				ja.push_back(j+1);
				continue;
			}
			if (i < j)
			{
				for (ind = matrix.m_rowptr[j]; ind < matrix.m_rowptr[j + 1]; ++ind)
					if (matrix.m_colind[ind] == i)
					{
						a.push_back(matrix.m_valU[ind]);
						ja.push_back(j+1);
						break;
					}
			}
			else
			{
				for (ind = matrix.m_rowptr[i]; ind < matrix.m_rowptr[i + 1]; ++ind)
					if (matrix.m_colind[ind] == j)
					{
						a.push_back(matrix.m_valL[ind]);
						ja.push_back(j+1);
						break;
					}
			}
		}
		ia.push_back(a.size()+1);
	}


	vector<int> perm(matrix.m_size);
	//perm[35] = 1;
	int nrhs = 1;
	int parm[64];
	double dparm[64];
	for (int i = 0; i < 64; ++i)
		parm[i] = 0;
	int msqlvl = 0;
	vector<double> b{ rhs };
	int error = 0;
	vector<double> sl(b.size());


	int size = matrix.m_size;
	double* _a = new double[a.size()];
	int* _ia = new int[ia.size()];
	int* _ja = new int[ja.size()];
	int* _perm = new int[perm.size()];
	double* _b = new double[b.size()];
	double* _sl = new double[sl.size()];
	for (int i = 0; i < a.size(); ++i)
		_a[i] = a[i];
	for (int i = 0; i < ia.size(); ++i)
		_ia[i] = ia[i];
	for (int j = 0; j < ja.size(); ++j)
		_ja[j] = ja[j];
	//for (int i = 0; i < perm.size(); ++i)
		//_perm[i] = perm[i];
	for (int i = 0; i < b.size(); ++i)
		_b[i] = b[i];
	for (int i = 0; i < sl.size(); ++i)
		_sl[i] = sl[i];
	int solver = 0;
	//pardisoinit(pt, &mtype, &solver);
	//pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, _a, _ia, _ja, _perm, &nrhs, parm, &msqlvl, _b, _sl, &error);
	phase = -1;
	//pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, _a, _ia, _ja, _perm, &nrhs, parm, &msqlvl, _b, _sl, &error);
	//cout << "error: " << error << endl;
	sol.resize(sl.size());
	for (int i = 0; i < sl.size(); ++i)
		sol[i] = _sl[i];
	delete[] _a;
	delete[] _ia;
	delete[] _ja;
	delete[] _perm;
	delete[] _b;
	delete[] _sl;
}
void ESolver::Gauss(Matrix& matrix, const std::vector<double>& rhs, std::vector<double>& sol)
{
	auto size = matrix.m_size;
	auto n = size;
	sol.resize(size);
	vector<vector<double>> A;
	A.resize(size);
	for (unsigned int i = 0; i < size; ++i)
	{
		A[i].resize(size + 1);
		for (unsigned int j = 0; j < size; ++j)
			A[i][j] = matrix.m_elem[i][j];
		A[i][size] = rhs[i];
	}
	for (int i = 0; i < size; ++i)
	{
		double max_elem = fabs(A[i][i]);
		int max_row = i;
		for (int j = i + 1; j < size; ++j)
			if (max_elem < fabs(A[j][i]))
			{
				max_elem = fabs(A[j][i]);
				max_row = j;
			}

		for (int j = i; j < size + 1; ++j)
		{
			double temp = A[max_row][j];
			A[max_row][j] = A[i][j];
			A[i][j] = temp;
		}

		for (int j = i + 1; j < size; ++j)
		{
			double temp = -A[j][i] / A[i][i];
			for (int k = i; k < size + 1; ++k)
				if (i == k)
					A[j][k] = 0;
				else
					A[j][k] += temp * A[i][k];
		}
	}
	for (int i = size - 1; i >= 0; --i)
	{
		sol[i] = A[i][size] / A[i][i];
		for (int j = i - 1; j >= 0; --j)
			A[j][size] -= A[j][i] * sol[i];
	}
}


void ESolver::Gauss(const Matrix& matrix, double* r)
{
	auto size = matrix.m_size;
	auto n = size;
	vector<vector<double>> A;
	A.resize(size);
	for (unsigned int i = 0; i < size; ++i)
	{
		A[i].resize(size + 1);
		for (unsigned int j = 0; j < size; ++j)
			A[i][j] = matrix.m_elem[i][j];
		A[i][size] = r[i];
	}
	for (int i = 0; i < size; ++i)
	{
		double max_elem = fabs(A[i][i]);
		int max_row = i;
		for (int j = i + 1; j < size; ++j)
			if (max_elem < fabs(A[j][i]))
			{
				max_elem = fabs(A[j][i]);
				max_row = j;
			}

		for (int j = i; j < size + 1; ++j)
		{
			double temp = A[max_row][j];
			A[max_row][j] = A[i][j];
			A[i][j] = temp;
		}

		for (int j = i + 1; j < size; ++j)
		{
			double temp = -A[j][i] / A[i][i];
			for (int k = i; k < size + 1; ++k)
				if (i == k)
					A[j][k] = 0;
				else
					A[j][k] += temp * A[i][k];
		}
	}
	for (int i = size - 1; i >= 0; --i)
	{
		r[i] = A[i][size] / A[i][i];
		for (int j = i - 1; j >= 0; --j)
			A[j][size] -= A[j][i] * r[i];
	}
}

void ESolver::Gauss(const Matrix& matrix, double* rhs, double* sol)
{
	auto size = matrix.m_size;
	auto n = size;
	vector<vector<double>> A;
	A.resize(size);
	for (unsigned int i = 0; i < size; ++i)
	{
		A[i].resize(size + 1);
		for (unsigned int j = 0; j < size; ++j)
			A[i][j] = matrix.m_elem[i][j];
		A[i][size] = rhs[i];
	}
	for (int i = 0; i < size; ++i)
	{
		double max_elem = fabs(A[i][i]);
		int max_row = i;
		for (int j = i + 1; j < size; ++j)
			if (max_elem < fabs(A[j][i]))
			{
				max_elem = fabs(A[j][i]);
				max_row = j;
			}

		for (int j = i; j < size + 1; ++j)
		{
			double temp = A[max_row][j];
			A[max_row][j] = A[i][j];
			A[i][j] = temp;
		}

		for (int j = i + 1; j < size; ++j)
		{
			double temp = -A[j][i] / A[i][i];
			for (int k = i; k < size + 1; ++k)
				if (i == k)
					A[j][k] = 0;
				else
					A[j][k] += temp * A[i][k];
		}
	}
	for (int i = size - 1; i >= 0; --i)
	{
		sol[i] = A[i][size] / A[i][i];
		for (int j = i - 1; j >= 0; --j)
			A[j][size] -= A[j][i] * sol[i];
	}
}

void ESolver::Gauss(const Matrix& matrix, const double* rhs, double* sol)
{
	auto size = matrix.m_size;
	auto n = size;
	vector<vector<double>> A;
	A.resize(size);
	for (unsigned int i = 0; i < size; ++i)
	{
		A[i].resize(size + 1);
		for (unsigned int j = 0; j < size; ++j)
			A[i][j] = matrix.m_elem[i][j];
		A[i][size] = rhs[i];
	}
	for (int i = 0; i < size; ++i)
	{
		double max_elem = fabs(A[i][i]);
		int max_row = i;
		for (int j = i + 1; j < size; ++j)
			if (max_elem < fabs(A[j][i]))
			{
				max_elem = fabs(A[j][i]);
				max_row = j;
			}

		for (int j = i; j < size + 1; ++j)
		{
			double temp = A[max_row][j];
			A[max_row][j] = A[i][j];
			A[i][j] = temp;
		}

		for (int j = i + 1; j < size; ++j)
		{
			double temp = -A[j][i] / A[i][i];
			for (int k = i; k < size + 1; ++k)
				if (i == k)
					A[j][k] = 0;
				else
					A[j][k] += temp * A[i][k];
		}
	}
	for (int i = size - 1; i >= 0; --i)
	{
		sol[i] = A[i][size] / A[i][i];
		for (int j = i - 1; j >= 0; --j)
			A[j][size] -= A[j][i] * sol[i];
	}
}

const std::vector<double> ESolver::Solve(MatrixDiag & matrix, const std::vector<double>& rhs, std::vector<double>& sol, std::vector<double>& residual, const int iter, const double eps)
{
	auto size{ matrix.GetSize() };
	m_solution.resize(size);
	// y = a * x + y
	// y = 0; a = 1./d; x = r
//#pragma omp parallel for shared(rhs, matrix.m_valDiag) private(i)
	for (int i = 0; i < size; ++i)
	{
		m_solution[i] = rhs[i] / matrix.m_valDiag[i];
	}
	return m_solution;
}
void ESolver::BiCGStabPrecond()
{
	LUPrec();
	double *r, *rt, *s, *p, *v, *t;
	double alp, ro, omega, omegal, rol;
	int k = 0;
	r = new double[m_matrix.m_size];
	rt = new double[m_matrix.m_size];
	p = new double[m_matrix.m_size];
	v = new double[m_matrix.m_size];
	s = new double[m_matrix.m_size];
	t = new double[m_matrix.m_size];
	ro = alp = omega = omegal = rol = 1;
	// v = p = 0;
	double *temp = new double[m_matrix.m_size];
	int i = 0;
	MatrixprodVector(temp, m_rightvector, m_matrix.m_valDiag, m_matrix.m_valL, m_matrix.m_valU, m_matrix.m_rowptr, m_matrix.m_colind, m_matrix.m_size);

	for (int i = 0; i < m_matrix.m_size; ++i)
	{
		r[i] = m_rightvector[i] - temp[i];
		p[i] = v[i] = 0;
	}
	LSolve(temp, r);
	for (int i = 0; i < m_matrix.m_size; ++i)
		rt[i] = r[i];
	double residual = Norm(r, m_matrix.m_size) / Norm(m_rightvector, m_matrix.m_size);
	while ((k < m_maxiter) && (residual > m_eps))
	{
		alp = DotProd(r, rt, m_matrix.m_size);
		//USolve()
		++k;
	}

}
void ESolver::LUPrec()
{
	m_LUvalL = new double[m_matrix.m_gsize];
	m_LUvalU = new double[m_matrix.m_gsize];
	m_LUvalDiag = new double[m_matrix.m_size];
	m_LUvalDiag[0] = m_matrix.m_valDiag[0];
	int i, j, k, l;
	for (i = 1; i < m_matrix.m_size; ++i)
	{
		m_LUvalDiag[i] = m_matrix.m_valDiag[i];
		for (j = m_matrix.m_rowptr[i]; j < m_matrix.m_rowptr[i + 1]; ++j)
		{
			m_LUvalU[j] = m_LUvalU[j];
			m_LUvalL[j] = m_matrix.m_valL[j];
			for (k = m_matrix.m_rowptr[i]; k < j; ++k)
				for (l = m_matrix.m_rowptr[m_matrix.m_colind[j]]; l < m_matrix.m_rowptr[m_matrix.m_colind[j]]; ++l)
				{
				if (m_matrix.m_colind[k] == m_matrix.m_colind[l])
				{
					m_LUvalU[j] -= m_LUvalU[k] * m_LUvalL[l];
					m_LUvalL[j] -= m_LUvalU[l] * m_LUvalL[k];
				}
				}
			if (fabs(m_LUvalDiag[m_matrix.m_colind[j]]) > m_eps)
				m_LUvalL[j] /= m_LUvalDiag[m_matrix.m_colind[j]];
			else
			{
				std::cout << "FAIL" << std::endl;
				return;
			}
			m_LUvalDiag[i] -= m_LUvalU[j] * m_LUvalL[j];
		}
	}
}
void ESolver::LSolve(double *f, double *res)
{
	double* temp = new double[m_matrix.m_size];
	for (int i = 0; i < m_matrix.m_size; ++i)
	{
		temp[i] = f[i];
		for (int j = m_matrix.m_rowptr[i]; j < m_matrix.m_rowptr[i + 1]; ++i)
		{
			temp[i] -= res[m_matrix.m_colind[j]] * m_LUvalL[j];
		}
		res[i] = temp[i] / m_LUvalDiag[i];
	}
	delete[] temp;
}
void ESolver::USolve(double *f, double *res)
{
	double* temp = new double[m_matrix.m_size];
	for (int i = 0; i < m_matrix.m_size; ++i)
		temp[i] = f[i];
	for (int i = m_matrix.m_size - 1; i > -1; --i)
	{
		res[i] = temp[i] / m_LUvalDiag[i];
		for (int j = m_matrix.m_rowptr[i]; j < m_matrix.m_rowptr[i + 1]; ++i)
		{
			temp[m_matrix.m_colind[j]] -= res[i] * m_LUvalU[j];
		}

	}
	delete[] temp;
}
ESolver::~ESolver()
{

}
