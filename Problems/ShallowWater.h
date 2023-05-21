#ifndef CORENC_PROBLEMS_SHALLOWWATER_H_
#define CORENC_PROBLEMS_SHALLOWWATER_H_

#include "Problems.h"
#include <vector>
#include "../CoreNCFEM/Parameter.h"
#include <map>
#include <tuple>
namespace corenc
{
	class CShallowWater : public CProblem
	{
		using boundary = std::tuple<int, Mesh::parameter<double>, Mesh::parameter<double>>;
	public:
		CShallowWater();
		~CShallowWater();
		Terms									getTerm(const unsigned int) const;
		const unsigned int						getNumberOfTerms() const;
		const int								setTerm(const unsigned int, const Terms&);
		const int								addTerm(const Terms&);
		const int								removeTerm(const Terms&);
		const int								load_parameters(const std::string& file_name);
		const double							get_parameter(const Terms&, const int element_type, const Mesh::Point&) const;
		const double							get_parameter(const Terms&, const int element_number, const int element_type, const Mesh::Point&) const;
		const double							get_boundary_parameter(const int type, const int element_type, const Mesh::Point&) const;
		const double							get_boundary_parameter(const int type, const int element_number, const int element_type, const Mesh::Point&) const;
		const int								get_number_of_boundaries() const;
		const double							get_solution(const int sys_number, const int element_type, const int element_number, const Mesh::Point&) const;
		const int								get_boundary_type(const int number) const;
		const int								add_parameter(const Terms&, const int element_type, const Mesh::parameter<double>& value);
		const int								set_parameter(const Terms&, const int element_type, const Mesh::parameter<double>& value);
		const int								set_boundary_parameter(const int type, const int element_type, const boundary& value);
		// 1st and 2nd types of boundary conditions
		const int								add_boundary_parameter(const int type, const int element_type, const Mesh::parameter<double>& value);
		// 3rd type of boundary conditions
		const int								add_boundary_parameter(const int element_type, const Mesh::parameter<double>& value, const Mesh::parameter<double>& value2);
	private:
		std::vector<Terms>						m_terms;
		std::map<int, Mesh::parameter<double>>	m_params;
		std::map<int, Mesh::parameter<double>>	m_srcs;
		std::map<int, boundary>					m_bounds;
		int										m_total_params;
		int										m_total_srcs;
		int										m_total_bounds;
	};
}
#endif // !CORENC_PROBLEMS_SHALLOWWATER_H_
