#include "DistributionFunction.hxx"

DistributionFunction::DistributionFunction(size_t nx, size_t ny, size_t nv) :
	m_nx(nx),
	m_ny(ny),
	m_nv(nv),
	ptr_f(new double[nx * ny * nv]) {}

DistributionFunction::~DistributionFunction() {}

double& DistributionFunction::operator()(size_t i, size_t j, size_t k)
{
	return ptr_f[(i * m_ny + j) * m_nv + k];
}