#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H
#include <vector>
#include <complex>
namespace corenc
{
    namespace solvers
    {
        template<class Matrix, class Solver>
        class eigen_solver
        {
        public:
            eigen_solver(){}
            ~eigen_solver(){}
            void        rayleigh(Matrix* A, Matrix* B, Solver* esl, std::complex<double>* mu0, double* x0, const int n) const
            {
                std::vector<std::complex<double>> x(n);
                std::vector<std::complex<double>> y(n);
                std::vector<std::complex<double>> lam(n);
                double norm_mu = 0;
                double norm_x = 0;
                for (int i = 0; i < n; ++i)
                {
                    norm_mu += std::norm(mu0[i]) * std::norm(mu0[i]);
                    norm_x += std::norm(x0[i]) * std::norm(x0[i]);
                }
                norm_mu = sqrt(norm_mu);
                norm_x = sqrt(norm_x);
                for (int i = 0; i < n; ++i)
                {
                    x[i] = x0[i] / norm_x;
                    y[i] = mu0[i] / norm_mu;
                }
                std::complex<double> temp(0, 0);
                temp =
            }
        };
    }
}
#endif // EIGEN_SOLVER_H
