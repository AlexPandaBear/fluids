#include "LBMKernel.hxx"

LBMKernel::LBMKernel(LBMConfig& config) :
	m_config(config),
{
	build_space();
	build_time();
	initialize_f();
}

LBMKernel::LBMKernel(std::string config_file) : LBMKernel(LBMConfig(config_file)) {}

LBMKernel::~LBMKernel() {}

void LBMKernel::build_space()
{
	m_space = std::vector<std::vector<double>>(0, std::vector<double>(0, 0.));

	for (size_t d = 0; d < m_config.get_nb_dim(); d++)
	{
		size_t n(m_config.get_grid(d));
		std::pair<double, double> bounds(m_config.get_space_bounds(d));

		double h((bounds.second - bounds.first) / (n-1));

		std::vector<double> v(n, 0.);
		
		for (size_t i = 0; i < n; i++)
		{
			v[i] = bounds.first + i*h;
		}

		m_space.push_back(v);
	}
}

void LBMKernel::build_time()
{
	std::pair<double, double> bounds(m_config.get_time_bounds());
	size_t n(m_config.get_nb_steps());

	m_time = std::vector<double>(n+1, 0.);

	double dt((bounds.second - bounds.start) / n);

	for (size_t t = 0; t < n+1; t++)
	{
		m_time[t] = bounds.first + t*dt;
	}
}

void LBMKernel::initialize_f()
{
	m_f = DistributionFunction(m_config.get_grid(), m_config.get_velocity_scheme());
}

void initialize_operators()
{
	m_collision = BGK(m_config.get_mu(), m_config.get_gamma());
	m_streaming = Stream2D(m_config.get_boundary_conditions());
}

void LBMKernel::set_initial_state(Field& rho, Field& U, Field& E)
{
	std::vector<size_t> sizes(m_f.get_sizes());

	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			m_f.set_equilibrium(i, j, rho({i,j}), U({i,j,0}), U({i,j,1}), E({i,j}), m_config.get_gamma());
		}
	}
}

void LBMKernel::simulate()
{
	m_streaming.initialize(m_f);

	for (size_t t = 0; t < nb_steps; t++)
	{
		m_collision(m_f, m_time[t+1] - m_time[t]);
		m_streaming(m_f);
	}
}