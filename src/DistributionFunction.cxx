#include "DistributionFunction.hxx"

DistributionFunction::DistributionFunction(size_t nx, size_t ny, size_t nv) :
	m_nx(nx),
	m_ny(ny),
	m_nv(nv),
	m_w(nv, 0.),
	m_ex(nv, 0.),
	m_ey(nv, 0.),
	ptr_f(new double[nx * ny * nv]) {}

DistributionFunction::DistributionFunction(size_t nx, size_t ny, size_t nv, std::string vs) :
	DistributionFunction(nx, ny, nv)
{
	if (vs == "D2Q9")																					//	4-----3-----2
	{																									//	| *   |   * |
		set_velocity_scheme({0, 1, 1, 0, -1, -1, -1, 0, 1},												//	|   * | *   |
							{0, 0, 1, 1, 1, 0, -1, -1, -1},												//	5-----0-----1
							{4./9., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36.});		//	|   * | *   |
	}																									//	| *   |   * |
																										//	6-----7-----8
	else
	{
		throw std::invalid_argument("Unknown velocity scheme");
	}
}

DistributionFunction::~DistributionFunction() {}

void DistributionFunction::set_velocity_scheme(std::vector<int> ex, std::vector<int> ey, std::vector<double> weights)
{
	for (size_t k = 0; k < m_nv; k++)
	{
		m_w[k] = weights[k];
		m_ex[k] = ex[k];
		m_ey[k] = ey[k];
		m_e2[k] = ex[k]*ex[k] + ey[k]*ey[k];
	}
}

void DistributionFunction::display_velocity_scheme() const
{
	std::cout << "Ex = [ ";

	for (size_t k = 0; k < m_nv; k++)
	{
		std::cout << m_ex[k] << " ";
	}

	std::cout << "]" << std::endl;

	std::cout << "Ey = [ ";

	for (size_t k = 0; k < m_nv; k++)
	{
		std::cout << m_ey[k] << " ";
	}

	std::cout << "]" << std::endl;

	std::cout << "We = [ ";

	for (size_t k = 0; k < m_nv; k++)
	{
		std::cout << m_w[k] << " ";
	}

	std::cout << "]" << std::endl;
}

size_t DistributionFunction::get_size_x() const
{
	return m_nx;
}

size_t DistributionFunction::get_size_y() const
{
	return m_ny;
}

size_t DistributionFunction::get_size_v() const
{
	return m_nv;
}

int DistributionFunction::get_ex(size_t k) const
{
	return m_ex[k];
}

int DistributionFunction::get_ey(size_t k) const
{
	return m_ey[k];
}

double DistributionFunction::get_w(size_t k) const
{
	return m_w[k];
}

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