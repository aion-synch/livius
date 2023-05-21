#include "Parameter.h"

using namespace std;
using namespace corenc;
using namespace Mesh;

CParameter::CParameter() :
	m_diffusion{ 0 },
	m_advection{ 0 },
	m_mass{ 0 }
{

}

//CParameter::CParameter(const double _diff, const double _adv, const double _mass) :
//	m_diffusion(_diff),
//	m_advection(_adv),
//	m_mass(_mass)
//{
//
//}

CParameter::CParameter(const parameter<double>& _diff, const parameter<double>& _adv, const parameter<double>& _mass) :
	m_diffusion(_diff),
	m_advection(_adv),
	m_mass(_mass)
{

}

CParameter::CParameter(const Parameters& type, const parameter<double>& p)
{
	switch (type)
	{
	case corenc::Parameters::DIFFUSION:
		m_diffusion = p;
		break;
	case corenc::Parameters::MASS:
		m_mass = p;
		break;
	case corenc::Parameters::ADVECTION:
		m_advection = p;
		break;
	default:
		break;
	}
}


CParameter::~CParameter()
{

}

const double CParameter::GetDiffusion() const
{
	return m_diffusion.get(Point(0,0,0));
}

const double CParameter::GetAdvection() const
{
	return m_advection.get(Point(0, 0, 0));
}

const double CParameter::GetMass() const
{
	return m_mass.get(Point(0, 0, 0));
}

const double CParameter::GetDiffusion(const Point& p) const
{
	return m_diffusion.get(p);
}

const double CParameter::GetAdvection(const Point& p) const
{
	return m_advection.get(p);
}

const double CParameter::GetMass(const Point& p) const
{
	return m_mass.get(p);
}