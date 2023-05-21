#pragma once
#ifndef CORENC_METHODS_FEANALYSIS_H_
#define CORENC_METHODS_FEANALYSIS_H_
#include <vector>
#include "../Point.h"
namespace corenc
{
	namespace method
	{
		template<class Method1, class Method2, class Mesh1, class Mesh2>
		class FEAnalysis
		{
		public:
			FEAnalysis() {};
			~FEAnalysis() {};
			const double					L2Norm(	const Method1& method1, 
													const Method2& method2, 
													const Mesh1& mesh1, 
													const Mesh2& mesh2, 
													const std::vector<double>& w1,
													const std::vector<double>& w2) const;
		};
		template<class Method1, class Method2, class Mesh1, class Mesh2>
		const double FEAnalysis<Method1, Method2, Mesh1, Mesh2>::L2Norm(
			const Method1& method1, 
			const Method2& method2, 
			const Mesh1& mesh1, 
			const Mesh2& mesh2, 
			const std::vector<double>& w1,
			const std::vector<double>& w2) const
		{
			double sum{ 0 }, sum2{0};
			double res, res2;
			int j;
			std::vector<int> dofs;
			int order = mesh1.GetElement(0)->GetDoFs();
			dofs.resize(order);
			std::vector<Mesh::Point> points(order);
			auto sub = [&](const Mesh::Point& p)
			{
				return (method1.GetValue(p, w1) - method2.GetValue(p, w2)) * (method1.GetValue(p, w1) - method2.GetValue(p, w2));
			};
			auto r = [&](const Mesh::Point& p)
			{
				return method1.GetValue(p, w1);
			};
			const int n = (int)mesh1.GetNumberOfElements();
			for (int i = 0; i < n; ++i)
			{
				const auto& elem = mesh1.GetElement(i);
				for (j = 0; j < order; ++j)
					points[j] = mesh1.GetNode(elem->GetNode(j));
				res = elem->Integrate(sub, points);
				res2 = elem->Integrate(r, points);
				sum += res;
				sum2 += res2;
			}
			if (dofs.size() > 0)
				std::vector<int>().swap(dofs);
			if (points.size() > 0)
				std::vector<Mesh::Point>().swap(points);
			return sqrt(sum/sum2);
		}
	}
}

#endif // !CORENC_METHODS_FEANALYSIS_H_

