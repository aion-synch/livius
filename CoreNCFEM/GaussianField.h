#ifndef CORENC_GAUSSIANPROCESS_H_
#define CORENC_GAUSSIANPROCESS_H_
#include <algorithm>
#include <vector>
#include <cmath>
#include "Point.h"
namespace corenc
{
	struct GaussianProcess
	{
		double sigma2;
		double l;
		double a;
		double b;
		double c;
		double A;
		double B;
		size_t K = 1;
		std::vector<double> lambda;
		GaussianProcess(const double L, const size_t num)
		{
			K = num;
			sigma2 = L;
			l = 2 * L;
			a = 1. / (4 * sigma2);
			b = 1. / (2 * l * l);
			c = sqrt(a * a + 2 * a * b);
			A = a + b + c;
			B = b / A;
			lambda.resize(K);
			for (size_t i = 0; i < K; ++i)
				lambda[i] = std::pow(B, i) * sqrt(2 * a / A);
		}
		const double He(const int i, const double x) const
		{
			switch (i)
			{
			case 0:
				return 1.;
			case 1:
				return x;
			case 2:
				return x * x - 1.;
			case 3:
				return x * x * x - 3. * x;
			case 4:
				return x * x * x * x - 6. * x * x + 3.;
			case 5:
				return x * x * x * x * x - 10. * x * x * x + 15. * x;
			case 6:
				return x * x * x * x * x * x - 15. * x * x * x * x + 45. * x * x - 25.;
			case 7:
				return x * x * x * x * x * x * x - 21. * x * x * x * x * x + 105. * x * x * x - 105. * x;
			case 8:
				return x * x * x * x * x * x * x * x - 28. * x * x * x * x * x * x + 210. * x * x * x * x - 420. * x * x + 105;
			case 9:
				return x * x * x * x * x * x * x * x * x - 36. * x * x * x * x * x * x * x + 378. * x * x * x * x * x - 1260. * x * x * x + 945. * x;
			default:
				return x * x * x * x * x * x * x * x * x * x - 45. * x * x * x * x * x * x * x * x + 630. * x * x * x * x * x * x - 3150. * x * x * x *x + 4725. * x * x - 945.;
				break;
			}
		};
		const double phi(const int i, const double x) const
		{
			return exp(-(c - a) * x * x) * He(i, x * sqrt(2 * c));
		};
	};
	/*enum class gkernels
	{
		gexponent,
		gker1,
		gker2,
		gker3
	};*/
	struct GaussianKernel
	{
		int N;
		const double gpexp(const Mesh::Point& a) const
		{
            return exp(-12.5 * (a.x * a.x + a.y * a.y));
		}

		const double gpstep(const Mesh::Point& a) const
		{
			if (fabs(a.x) < 0.5 && fabs(a.y) < 0.5)
				return 1.;
			return 0.;
		}
		std::vector<Mesh::Point> _centrs;
		GaussianKernel(const int _n, const std::vector<Mesh::Point>& centers) :
			N{ _n }, _centrs{ centers } {}
		const double get_gp(const std::vector<double>& a, const Mesh::Point& p) const
		{
			double sum = 0;
			for (auto i = 0; i < N; ++i)
				sum += a[i] * gpexp(p - _centrs[i]);
			return sum;
		}

	};
}
#endif // CORENC_GAUSSIANPROCESS_H_
