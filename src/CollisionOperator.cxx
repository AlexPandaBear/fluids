#include "CollisionOperator.hxx"

void CollisionOperator::set_f(DistributionFunction const& f)
{
	m_f = f;
}

BGK::BGK(double nu, double r, std::function<double(..)> f_eq) :
	m_nu(nu),
	m_r(r),
	m_f_eq(f_eq) {}

BGK::~BGK() {}

double BGK::operator()(size_t i, size_t j, size_t k) const
{
	return 0.;
}