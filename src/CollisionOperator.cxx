#include "CollisionOperator.hxx"

BGK::BGK(double nu, double gamma, double dt) :
	m_nu(nu),
	m_gamma(gamma),
	m_dt(dt) {}

BGK::~BGK() {}

std::vector<double> BGK::operator()(DistributionFunction const& f, size_t i, size_t j) const
{
	double rho(f.compute_rho(i, j));
	
	double ux(f.compute_rho_ux(i, j) / rho);
	double uy(f.compute_rho_uy(i, j) / rho);

	double c2(f.compute_rho_e(i, j) * (m_gamma-1.) / rho);
	double inv_c2(1. / c2);

	double u2_c2((ux*ux + uy*uy) * inv_c2);
	double inv_tau = c2 / (m_nu + 0.5 * m_dt * c2);

	std::vector<double> res(m_nv, 0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		double eu_c2((m_ex[k]*ux + m_ey[k]*uy) * inv_c2);
		res[k] = inv_tau * (m_w[k] * rho * (1. + eu_c2 + 0.5*(eu_c2*eu_c2 - u2_c2)) - f(i, j, k));
	}

	return res;
}