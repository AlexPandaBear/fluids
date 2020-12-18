#include "CollisionOperator.hxx"

CollisionOperator::~CollisionOperator() {}

void CollisionOperator::operator()(DistributionFunction& f, double dt) const {}

BGK::BGK(double nu, double gamma) :
	m_nu(nu),
	m_gamma(gamma) {}

BGK::~BGK() {}

void BGK::operator()(DistributionFunction& f, double dt) const
{
	size_t nx(f.get_size_x());
	size_t ny(f.get_size_y());
	size_t nv(f.get_size_v());

	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			double rho(f.compute_rho(i, j));
			
			double ux(f.compute_rho_ux(i, j) / rho);
			double uy(f.compute_rho_uy(i, j) / rho);

			double c2(f.compute_rho_e(i, j) * (m_gamma-1.) / rho);
			double inv_c2(1. / c2);

			double u2_c2((ux*ux + uy*uy) * inv_c2);
			double inv_tau = c2 / (m_nu + 0.5 * dt * c2);


			for (size_t k = 0; k < nv; k++)
			{
				double eu_c2((f.get_ex(k) * ux + f.get_ey(k) * uy) * inv_c2);
				f(i, j, k) += dt * inv_tau * (f.get_w(k) * rho * (1. + eu_c2 + 0.5*(eu_c2*eu_c2 - u2_c2)) - f(i, j, k));
			}
		}
	}
}