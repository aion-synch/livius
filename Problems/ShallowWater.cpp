#include "ShallowWater.h"
#include <vector>
#include <istream>
#include <iostream>
#include <fstream>
using namespace corenc;

CShallowWater::CShallowWater()
{
	m_terms.resize(1);
	m_terms[0] = Terms::IDUDV;
	m_total_bounds = 0;
	m_total_params = 0;
	m_total_srcs = 0;
}

const int CShallowWater::removeTerm(const Terms& term)
{
	for (unsigned i{ 0 }; i < m_terms.size(); ++i)
		if (m_terms[i] == term)
			m_terms.erase(m_terms.begin() + i);
	return 0;
}

Terms CShallowWater::getTerm(const unsigned int i) const
{
	if (i < m_terms.size())
		return m_terms[i];
	return m_terms[0];
}

const unsigned int CShallowWater::getNumberOfTerms() const
{
	return m_terms.size();
}

const int CShallowWater::setTerm(const unsigned int i, const Terms& term)
{
	if (i < m_terms.size())
	{
		m_terms[i] = term;
		return 0;
	}
	return 1;
}

const int CShallowWater::addTerm(const Terms& term)
{
	m_terms.push_back(term);
	return 0;
}

const int CShallowWater::load_parameters(const std::string& file_name)
{
	std::ifstream ifs;
	ifs.open(file_name);
	if (!ifs.is_open())
		return 1;
	double diff, temp;
	int type, temp2, n;
	// material properties
	// int - number of materials
	// int, double - material, its conductivity
	// boundary conditions
	// int - number of conditions
	// int_2, int_1, double_1,.. - if int_1 < 3 -> 1st or 2nd type of conditions
	// if int_1 == 3 -> 3rd type of condtions
	// int_2 - type of boundaries
	// double_1 boundary condition
	// double_2 only for the 3rd type
	ifs >> n;
	if (n < 1)
		return 1;
	m_total_params = n;
	for (int i = 0; i < n; ++i)
	{
		ifs >> type >> diff;
		m_params[type] = Mesh::parameter<double>(diff);
	}
	ifs >> n;
	if (n < 0)
		return 1;
	m_total_bounds = n;
	for (int i = 0; i < n; ++i)
	{
		ifs >> type >> temp2;
		if (temp2 < 3)
		{
			ifs >> temp;
			m_bounds[type] = std::make_tuple(temp2, Mesh::parameter<double>(temp), Mesh::parameter<double>(0.));
		}
		else if (temp2 == 3)
		{
			ifs >> temp >> diff;
			m_bounds[type] = std::make_tuple(temp2, Mesh::parameter<double>(temp), Mesh::parameter<double>(diff));
		}
		else
			return 1;
	}
	return 0;
}

const double CShallowWater::get_parameter(const Terms& term, const int type, const Mesh::Point& p) const
{
	if (m_params.find(type) == m_params.end())
		return 0.;
	if (term == Terms::IDUDV)
		return m_params.find(type)->second.get(p);
	if (term == Terms::EFV)
		return m_params.find(type)->second.get(p);
	return 0.;
}

const double CShallowWater::get_parameter(const Terms& term, const int number, const int type, const Mesh::Point& p) const
{
	if (m_params.find(type) == m_params.end())
		return 0.;
	if (term == Terms::IDUDV)
		return m_params.find(type)->second.get(number, p);
	if (term == Terms::EFV)
		//return m_params.find(type)->second.get(number, p);
		return m_srcs.find(type)->second.get(number, p);
	return 0.;
}

const double CShallowWater::get_boundary_parameter(const int type, const int element_type, const Mesh::Point & p) const
{
	if (m_bounds.find(element_type) == m_bounds.end())
		return 0.0;
	if (type == 0)
		return std::get<1>(m_bounds.find(element_type)->second).get(p);
	return std::get<2>(m_bounds.find(element_type)->second).get(p);
}

const double CShallowWater::get_boundary_parameter(const int type, const int element_number, const int element_type, const Mesh::Point & p) const
{
	if (m_bounds.find(element_type) == m_bounds.end())
		return 0.0;
	if (type == 0)
		return std::get<1>(m_bounds.find(element_type)->second).get(element_number, p);
	return std::get<2>(m_bounds.find(element_type)->second).get(element_number, p);
}

const int CShallowWater::get_number_of_boundaries() const
{
	return m_total_bounds;
}

const double CShallowWater::get_solution(const int sys_number, const int element_type, const int element_number, const Mesh::Point &) const
{
	return 0.0;
}

const int CShallowWater::get_boundary_type(const int number) const
{
	auto it = m_bounds.begin();
	int i = 0;
	for (; it != m_bounds.end() && i < number; ++it, ++i);
	return it->first;
}

const int CShallowWater::add_parameter(const Terms & term, const int element_type, const Mesh::parameter<double>& value)
{
	switch (term)
	{
	case Terms::IDUDV:
		m_params[element_type] = value;
		++m_total_params;
		return 0;
	case Terms::EFV:
		m_srcs[element_type] = value;
		++m_total_srcs;
		return 0;
	default:
		break;
	}
	return 1;
}

const int CShallowWater::set_parameter(const Terms & term, const int element_type, const Mesh::parameter<double>& value)
{
	switch (term)
	{
	case Terms::IDUDV:
		m_params[element_type] = value;
		return 0;
	default:
		break;
	}
	return 1;
}

const int CShallowWater::set_boundary_parameter(const int type, const int element_type, const boundary & value)
{
	m_bounds[element_type] = value;
	return 0;
}

const int CShallowWater::add_boundary_parameter(const int type, const int element_type, const Mesh::parameter<double> & value)
{
	m_bounds[element_type] = std::make_tuple(type, value, 0);
	++m_total_bounds;
	return 0;
}

const int CShallowWater::add_boundary_parameter(const int element_type, const Mesh::parameter<double> & value, const Mesh::parameter<double> & value2)
{
	m_bounds[element_type] = std::make_tuple(3, value, value2);
	++m_total_bounds;
	return 0;
}

CShallowWater::~CShallowWater()
{
	if (m_terms.size() > 0)
		std::vector<Terms>().swap(m_terms);
	if (m_bounds.size() > 0)
		std::map<int, boundary>().swap(m_bounds);
	if (m_params.size() > 0)
		std::map<int, Mesh::parameter<double>>().swap(m_params);
}
