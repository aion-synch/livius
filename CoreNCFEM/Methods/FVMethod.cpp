#include "FVMethod.h"
using namespace std;
using namespace corenc;
using namespace method;
using namespace Mesh;


FVMethod1d::FVMethod1d()
{

}

FVMethod1d::~FVMethod1d()
{

}

const int FVMethod1d::Solve(CMesh<CFESolution>* mesh, const std::function<const double(const double)>& flux_func, const FVFlux& flux_type, std::vector<double>& solution, const double time_step)
{
	auto size = mesh->getSolution().size();
	std::vector<double> sol(size);
	auto numerical_flux = [&](const double ul, const double ur)
	{
		return (flux_func(ul) + flux_func(ur)) / 2 - fmax(ul, ur) / 2 * (ur - ul);
	};
	auto sz = size - 2;
	auto min_size = mesh->getMinSize();
	for (int i = 0; i < sz; ++i)
		sol[i + 1] = mesh->getSolution()[i + 1] + time_step / min_size * (numerical_flux(mesh->getSolution()[i], mesh->getSolution()[i + 1]) - numerical_flux(mesh->getSolution()[i + 1], mesh->getSolution()[i + 2]));
	sol[size - 1] = sol[size - 2];
	sol[0] = sol[1];
	solution.resize(sol.size());
	solution = sol;
	mesh->updateSolution(sol);
	return 0;
}

const double FVMethod1d::GetSolution(const CMesh1D& g, const Point& p)
{
	return g.getSolution(g.FindElement(p), 0);
}