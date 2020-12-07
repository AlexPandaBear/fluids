#include "DistributionFunction.hxx"

DistributionFunction::DistributionFunction(size_t nx, size_t ny, size_t nv) :
	m_nx(nx),
	m_ny(ny),
	m_nv(nv),
	m_w(nv, 0.),
	m_ex(nv, 0.),
	m_ey(nv, 0.),
	ptr_f(new double[nx * ny * nv]) {}

DistributionFunction::~DistributionFunction() {}

double DistributionFunction::compute_rho(size_t i, size_t j) const
{
	double rho(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		rho += m_w[k] * ptr_f[(i * m_ny + j) * m_nv + k];
	}

	return rho;
}

double DistributionFunction::compute_rho_ux(size_t i, size_t j) const
{
	double ux(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		ux += m_w[k] * m_ex[k] * ptr_f[(i * m_ny + j) * m_nv + k];
	}

	return ux;
}

double DistributionFunction::compute_rho_uy(size_t i, size_t j) const
{
	double uy(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		uy += m_w[k] * m_ey[k] * ptr_f[(i * m_ny + j) * m_nv + k];
	}

	return uy;
}

double DistributionFunction::compute_rho_e(size_t i, size_t j) const
{
	double res(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		res += m_w[k] * m_e2[k] * ptr_f[(i * m_ny + j) * m_nv + k];
	}

	return 0.5*res;
}

double& DistributionFunction::operator()(size_t i, size_t j, size_t k)
{
	return ptr_f[(i * m_ny + j) * m_nv + k];
}