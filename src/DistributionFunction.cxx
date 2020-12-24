#include "DistributionFunction.hxx"

DistributionFunction::DistributionFunction(std::vector<size_t>, std::string vs) :
	m_nb_dim(sizes.size()),
	m_sizes(sizes),
	m_nv(0),
	m_e(0, std::vector<int>(0, 0)),
	m_e2(0, 0.),
	m_w(0, 0.),
	ptr_f(NULL)
{
	switch (vs)
	{
		case "D2Q9":
			set_velocity_scheme({{0, 0},						//	4-----3-----2
								 {1, 0},						//	| *   |   * |
								 {1, 1},						//	|   * | *   |
								 {0, 1},						//	5-----0-----1
								 {-1, 1},						//	|   * | *   |
								 {-1, 0},						//	| *   |   * |
								 {-1, -1},						//	6-----7-----8
								 {0, -1},
								 {1, -1}},
								{4./9., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36., 1./9., 1./36.});
		default:
			throw std::invalid_argument("Unknown velocity scheme");
	}
}

DistributionFunction::~DistributionFunction() {}

void DistributionFunction::set_velocity_scheme(std::vector<std::vector<int>> elem_vel, std::vector<double> weights)
{
	m_nv = weights.size();
	m_e = elem_vel;
	m_e2 = std::vector(m_nv, 0.)

	for (size_t k = 0; k < m_nv; k++)
	{
		for (size_t d = 0; d < m_nb_dim; d++)
		{
			m_e2[k] += pow(m_e[k][d], 2);
		}
	}

	m_w = weights;

	size_t size_tot(m_nv);
	for (size_t s : m_sizes)
	{
		size_tot *= s;
	}

	ptr_f.reset(new double[size_tot]);
}
/*
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
*/
std::vector<size_t> DistributionFunction::get_sizes() const
{
	return m_sizes;
}

std::vector<int> DistributionFunction::get_e(size_t k) const
{
	return m_e[k];
}

double DistributionFunction::get_w(size_t k) const
{
	return m_w[k];
}

std::vector<double> DistributionFunction::get_f(std::vector<size_t> pt) const
{
	size_t index(pt[0]);

	for (size_t d = 1; d < m_nb_dim; d++)
	{
		index *= m_sizes[i];
		index += pt[i];
	}

	index *= m_nv;

	std::vector<double> res(m_nv, 0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		res[k] = ptr_f[index + k];
	}

	return res;
}

double DistributionFunction::compute_rho(std::vector<size_t> pt) const
{
	std::vector<double> f(get_f(pt));
	double rho(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		rho += m_w[k] * f[k];
	}

	return rho;
}

double DistributionFunction::compute_rho_ui(std::vector<size_t> pt, size_t axis) const
{
	std::vector<double> f(get_f(pt));
	double res(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		res += m_w[k] * m_e[k][axis] * f[k];
	}

	return res;
}

double DistributionFunction::compute_rho_e(std::vector<size_t> pt) const
{
	std::vector<double> f(get_f(pt));
	double res(0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		res += m_w[k] * m_e2[k] * f[k];
	}

	return 0.5*res;
}

std::vector<double> DistributionFunction::compute_macro_variables(std::vector<size_t> pt) const
{
	std::vector<double> f(get_f(pt));
	std::vector<double> res(m_nb_dim + 2, 0.);

	for (size_t k = 0; k < m_nv; k++)
	{
		res[0] += m_w[k] * f[k];
	}

	for (size_t ax = 0; ax < m_nb_dim; ax++)
	{
		for (size_t k = 0; k < m_nv; k++)
		{
			res[1 + ax] += m_w[k] * m_e[k][ax] * f[k];
		}

		res[1 + ax] /= res[0];
	}

	for (size_t k = 0; k < m_nv; k++)
	{
		res[m_nb_dim + 1] += m_w[k] * m_e2[k] * f[k];
	}

	res[m_nb_dim + 1] /= res[0];

	return res;
}

void DistributionFunction::set_equilibrium(std::vector<size_t> pt, double rho, std::vector<double> u, double e, double gamma)
{
	double c2(e * (gamma-1.));
	double inv_c2(1. / c2);

	double u2(0.);

	for (size_t d = 0; d < m_nb_dim; d++)
	{
		u2 += pow(u[d], 2);
	}

	double u2_c2(u2 * inv_c2);
	double eu_c2;

	for (size_t k = 0; k < m_nv; k++)
	{
		eu_c2 = 0.;

		for (size_t d = 0; d < m_nb_dim; d++)
		{
			eu_c2 += m_e[k][d] * u[d];
		}

		eu_c2 *= inv_c2;

		(*this)(pt, k) += get_w(k) * rho * (1. + eu_c2 + 0.5*(eu_c2*eu_c2 - u2_c2));
	}
}

double& DistributionFunction::operator()(std::vector<size_t> pt, size_t k)
{
	size_t index(pt[0]);

	for (size_t d = 1; d < m_nb_dim; d++)
	{
		index *= m_sizes[i];
		index += pt[i];
	}

	index *= m_nv;

	return ptr_f[index + k];
}