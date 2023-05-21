#ifndef CORENC_METHODS_FINITEVOLUME_H_
#define CORENC_METHODS_FINITEVOLUME_H_

#include "../Grids/Mesh1D.h"

namespace corenc
{
	namespace method
	{
		enum class FVFlux
		{
			LaxFriedrichs,
			Upwind,
			Central,
			NOFLUX,
		};
		class FVMethod1d
		{
		public:
			FVMethod1d();
			~FVMethod1d();
			static const int						Solve(Mesh::CMesh<CFESolution>* mesh,
													const std::function<const double(const double)>& flux_func, 
													const FVFlux& flux_type, 
													std::vector<double>& new_solution, 
													const double time_step);
			static const double	 					GetSolution(const Mesh::CMesh1D& g, const Mesh::Point& p);
		};
	}
}
#endif // CORENC_METHODS_FINITEVOLUME_H_
