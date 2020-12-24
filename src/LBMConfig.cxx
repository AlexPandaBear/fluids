#include "LBMConfig.hxx"

LBMConfig::LBMConfig() :
	m_sim_name("undefined"),
	m_sim_type("undefined"),
	m_turb_model("undefined"),
	m_mu(0.),
	m_lambda(0.),
	m_r(0.),
	m_gamma(0.),
	m_nb_dim(0),
	m_space_bounds(0),
	m_space_grid(0),
	m_velocity_scheme("undefined"),
	m_time_bounds(0., 0.),
	m_nb_steps(0),
	m_output_dt(0.),
	m_output_type("undefined") {}

LBMConfig::LBMConfig(std::string config_file);

LBMConfig::~LBMConfig() {}

void LBMConfig::set_simulation_name(std::string name)
{
	m_sim_name = name;
}

void LBMConfig::set_simulation_type(std::string type)
{
	m_sim_type = type;
}

void LBMConfig::set_turbulence_model(std::string model)
{
	m_turb_model = model;
}

void LBMConfig::set_fluid_properties(double mu, double lambda, double r, double gamma)
{
	m_mu = mu;
	m_lambda = lambda;
	m_r = r;
	m_gamma = gamma;
}

void LBMConfig::set_space(std::vector<std::pair<double, double>> bounds, std::vector<size_t> grid)
{
	if (bounds.size() != grid.size())
	{
		throw std::invalid_argument("Arguments' sizes do not match");
	}

	m_nb_dim = bounds.size();
	m_space_bounds = bounds;
	m_space_grid = grid;
}

void LBMConfig::set_velocity_scheme(std::string scheme)
{
	m_velocity_scheme = scheme;
}

void LBMConfig::set_time(std::pair<double, double> bounds, size_t nb_steps)
{
	m_time_bounds = bounds;
	m_nb_steps = nb_steps;
}

void LBMConfig::set_output_parameters(double dt, std::string type)
{
	m_output_dt = dt;
	m_output_type = type;
}

std::string LBMConfig::get_simulation_name() const
{
	return m_sim_name;
}

std::string LBMConfig::get_simulation_type() const
{
	return m_sim_type;
}

std::string LBMConfig::get_turbulence_model() const
{
	return m_turb_model;
}

double get_mu() const
{
	return m_mu;
}

double get_lambda() const
{
	return m_lambda;
}

double LBMConfig::get_r() const
{
	return m_r;
}

double LBMConfig::get_gamma() const
{
	return m_gamma;
}

size_t LBMConfig::get_nb_dim() const
{
	return m_nb_dim;
}

std::pair<double, double> LBMConfig::get_space_bounds(size_t dimension) const
{
	return m_space_bounds[dimension];
}

size_t LBMConfig::get_grid(size_t dimension) const
{
	return m_space_grid[dimension];
}

std::vector<size_t> LBMConfig::get_grid() const
{
	return m_space_grid;
}

std::string LBMConfig::get_velocity_scheme() const
{
	return m_velocity_scheme;
}

std::pair<double, double> LBMConfig::get_time_bounds() const
{
	return m_time_bounds;
}

size_t LBMConfig::get_nb_steps() const
{
	return m_nb_steps;
}

double LBMConfig::get_output_dt() const
{
	return m_output_dt;
}

std::string LBMConfig::get_output_type() const
{
	return m_output_type;
}

void LBMConfig::save(std::string file_name) const;