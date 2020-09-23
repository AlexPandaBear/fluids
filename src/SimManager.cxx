#include "SimManager.hxx"

SimManager::SimManager() {}
SimManager::~SimManager() {}

void SimManager::defineTimeParameters(double tMax, size_t nb_steps)
{
	m_tMax = tMax;
	m_nb_steps = nb_steps;
}

void SimManager::defineGridParameters(double Lx, double Ly, size_t nx, size_t ny)
{
	m_Lx = Lx;
	m_Ly = Ly;
	m_nx = nx;
	m_ny = ny;
}

void SimManager::defineFluidProperties(double lambda, double rho, double cv, double mu)
{
	m_lambda = lambda;
	m_rho = rho;
	m_cv = cv;
	m_mu = mu;
}

void SimManager::defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> P0, std::vector<std::vector<double>> T0)
{
	m_data.reset_size(m_nb_steps, m_nx, m_ny);

	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_data.setTemperatureAt(0, i, j, T0[i][j]);
			m_data.setPressureAt(0, i, j, P0[i][j]);
			m_data.setXVelocityAt(0, i, j, U0[i][j]);
			m_data.setYVelocityAt(0, i, j, V0[i][j]);
		}
	}
}

void SimManager::defineBoundaryConditions(std::string type, double value)
{
	if (type == "temp")
	{
		m_temp_BC = true;
		m_flux_BC = false;
		m_BC_value = value;
	}
	else if (type == "flux")
	{
		m_temp_BC = false;
		m_flux_BC = true;
		m_BC_value = value;
	}
	else
	{
		throw std::invalid_argument("Bad boundary condition type");
	}
}

void SimManager::defineErrorTolerance(double press_err_max)
{
	m_press_err_max = press_err_max;
}


void SimManager::launchSimulation()
{
	IncompressibleKernel flow_kernel(m_data, m_tMax, m_nb_steps, m_Lx, m_Ly, m_nx, m_ny, m_rho, m_mu, m_press_err_max);
	flow_kernel.simulate();
	ThermalKernel th_kernel(m_data, m_tMax, m_nb_steps, m_Lx, m_Ly, m_nx, m_ny, m_lambda, m_rho, m_cv, m_temp_BC, m_flux_BC, m_BC_value);
	th_kernel.simulate();
}

double SimManager::getTemperatureAt(size_t t, size_t i, size_t j) const
{
	return m_data.getTemperatureAt(t, i, j);
}

double SimManager::getPressureAt(size_t t, size_t i, size_t j) const
{
	return m_data.getPressureAt(t, i, j);
}

double SimManager::getXVelocityAt(size_t t, size_t i, size_t j) const
{
	return m_data.getXVelocityAt(t, i, j);
}

double SimManager::getYVelocityAt(size_t t, size_t i, size_t j) const
{
	return m_data.getYVelocityAt(t, i, j);
}

void SimManager::setTemperatureAt(size_t t, size_t i, size_t j, double T)
{
	m_data.setTemperatureAt(t, i, j, T);
}

void SimManager::setPressureAt(size_t t, size_t i, size_t j, double P)
{
	m_data.setPressureAt(t, i, j, P);
}

void SimManager::setXVelocityAt(size_t t, size_t i, size_t j, double U)
{
	m_data.setXVelocityAt(t, i, j, U);
}

void SimManager::setYVelocityAt(size_t t, size_t i, size_t j, double V)
{
	m_data.setYVelocityAt(t, i, j, V);
}

void SimManager::saveData(std::string file_name) const
{
	m_data.saveData(file_name);
}

void SimManager::loadData(std::string file_name)
{
	m_data.loadData(file_name);
}

/*
double SimManager::get_tMax() const
{
	return m_tMax;
}

size_t SimManager::get_nb_steps() const
{
	return m_nb_steps;
}

double SimManager::get_dt() const
{
	return m_dt;
}

double SimManager::get_Lx() const
{
	return m_Lx;
}

double SimManager::get_Ly() const
{
	return m_Ly;
}

size_t SimManager::get_nx() const
{
	return m_nx;
}

size_t SimManager::get_ny() const
{
	return m_ny;
}

double SimManager::get_dx() const
{
	return m_dx;
}

double SimManager::get_dy() const
{
	return m_dy;
}

double SimManager::get_lambda() const
{
	return m_lambda;
}

double SimManager::get_rho() const
{
	return m_rho;
}

double SimManager::get_cv() const
{
	return m_cv;
}

double SimManager::get_mu() const
{
	return m_mu;
}

double SimManager::get_nu() const
{
	return m_nu;
}

double SimManager::get_Cv() const
{
	return m_Cv;
}

bool SimManager::get_temp_BC() const
{
	return m_temp_BC;
}

bool SimManager::get_flux_BC() const
{
	return m_flux_BC;
}

double SimManager::get_BC_value() const
{
	return m_BC_value;
}

double SimManager::get_press_err_max() const
{
	return m_press_err_max;
}
*/