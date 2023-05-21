#ifndef CORENC_PROBLEMS_DIFFUSIONSCALAR_H_
#define CORENC_PROBLEMS_DIFFUSIONSCALAR_H_

#include "Problems.h"
#include <vector>
#include "../CoreNCFEM/Parameter.h"
#include <map>
#include <tuple>
namespace corenc
{
	class CDiffusionScalar : public CProblem
	{
		using boundary = std::tuple<int, Mesh::parameter<double>, Mesh::parameter<double>>;
	public:
		CDiffusionScalar();
		~CDiffusionScalar();
		Terms									getTerm(const unsigned int) const;
		const unsigned int						getNumberOfTerms() const;
		const int								findTerm(const Terms&) const;
		const int								setTerm(const unsigned int, const Terms&);
		const int								addTerm(const Terms&);
		const int								removeTerm(const Terms&);
		const int								load_parameters(const std::string& file_name);
		const double							get_parameter(const Terms&, const int element_type, const Mesh::Point&) const;
		const double							get_parameter(const Terms&, const int element_number, const int element_type, const Mesh::Point&) const;
		const Mesh::Point						get_parameter(const Terms&, const int element_number, const int element_type, const Mesh::Point&, const int) const;
		const double							get_parameter(const Terms&, const int element_type, const int element_number, const int node, const Mesh::Point&) const;
		const Mesh::Point						get_parameter(const Terms&, const int element_type, const int element_number, const int node, const Mesh::Point&, const int v) const;
		const double							get_boundary_parameter(const int type, const int element_type, const Mesh::Point&) const;
		const double							get_boundary_parameter(const int type, const int element_type, const int element_number, const Mesh::Point&) const;
		const double							get_boundary_parameter(const int type, const int element_type, const int element_number, const int node, const Mesh::Point&) const;
		const int								get_number_of_boundaries() const;
		const int								get_boundary_type(const int number) const;
		const int								add_parameter(const Terms&, const int element_type, const double& value);
		const int								add_parameter(const Terms&, const int element_type, const Mesh::parameter<double>& value);
		const int								add_parameter(const Terms&, const int element_type, const Mesh::parameter<Mesh::Point>& value);
		const int								set_parameter(const Terms&, const int element_type, const Mesh::parameter<double>& value);
		const int								set_parameter(const Terms&, const int element_type, const Mesh::parameter<Mesh::Point>& value);
		const int								set_boundary_parameter(const int type, const int element_type, const boundary& value);
												// 1st and 2nd types of boundary conditions
		const int								add_boundary_parameter(const int type, const int element_type, const Mesh::parameter<double>& value);
												// 3rd type of boundary conditions
		const int								add_boundary_parameter(const int element_type, const Mesh::parameter<double>& value, const Mesh::parameter<double>& value2);
        const Mesh::point_source<double>        get_point_source(const int number) const;
        void                                    set_point_source(const int number, const Mesh::point_source<double>&);
        const int                               get_total_sources() const;
	private:
        std::vector<Terms>                          m_terms;
        std::map<int, Mesh::parameter<double>>      m_params;
		std::map<int, Mesh::parameter<Mesh::Point>>	m_vels;
        std::map<int, Mesh::parameter<double>>      m_srcs;
        std::map<int, Mesh::parameter<double>>      m_gams;
        std::map<int, boundary>                     m_bounds;
        //std::map<int, Mesh::point_source<double>>   m_pointsrcs;
        std::vector<Mesh::point_source<double>>     m_pointsrcs;
        int                                         m_total_params;
        int                                         m_total_srcs;
        int                                         m_total_gams;
        int                                         m_total_bounds;
	};
}
#endif // !CORENC_PROBLEMS_DIFFUSIONSCALAR_H_
