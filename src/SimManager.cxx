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

void SimManager::defineFluidProperties(double lambda, double rho, double cp, double mu)
{
	m_lambda = lambda;
	m_rho = rho;
	m_cp = cp;
	m_mu = mu;
}

void SimManager::defineBody(std::vector<std::vector<bool>> body)
{
	if (body.size() != m_nx || body[0].size() != m_ny)
	{
		throw std::invalid_argument("Body grid incompatible with meshing");
	}

	m_body = body;
}

void SimManager::defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> T0)
{
	m_data.reset_size(m_nb_steps, m_nx, m_ny);

	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_data.setTemperatureAt(0, i, j, T0[i][j]);
			m_data.setXVelocityAt(0, i, j, U0[i][j]);
			m_data.setYVelocityAt(0, i, j, V0[i][j]);
		}
	}
}

void SimManager::defineUniformDynamicBoundaryConditions(double U_bc, double V_bc, double P_bc)
{
	m_U_BC = std::vector<double>(2*m_nx*m_ny, U_bc);
	m_V_BC = std::vector<double>(2*m_nx*m_ny, V_bc);
	m_P_BC = std::vector<double>(2*m_nx*m_ny, P_bc);
}

void SimManager::defineDynamicBoundaryConditions(std::vector<double> v_U_bc, std::vector<double> v_V_bc, std::vector<double> v_P_bc)
{
	m_U_BC = v_U_bc;
	m_V_BC = v_V_bc;
	m_P_BC = v_P_bc;
}

void SimManager::defineThermalBoundaryConditions(std::string type, double value)
{
	if (type == "temp")
	{
		m_temp_BC = true;
		m_flux_BC = false;
		m_T_BC = value;
	}
	else
	{
		throw std::invalid_argument("Bad boundary condition type");
	}
}

void SimManager::defineFlowIntegrationParameters(double theta, double accuracy_u, double accuracy_v, double accuracy_p, bool clean_pressure)
{
	if (theta < 0. || theta > 1.)
	{
		throw std::invalid_argument("Theta parameter must be in [0,1]");
	}
	
	if (accuracy_u <= 0. || accuracy_v <= 0. || accuracy_p <= 0.)
	{
		throw std::invalid_argument("Accuracy parameters must be greater than 0");
	}
	
	m_theta_inc = theta;
	m_accuracy_u = accuracy_u;
	m_accuracy_v = accuracy_v;
	m_accuracy_p = accuracy_p;
	m_clean_pressure = clean_pressure;
}

void SimManager::defineThermalIntegrationParameters(double theta, double accuracy)
{
	if (theta < 0. || theta > 1.)
	{
		throw std::invalid_argument("Theta parameter must be in [0,1]");
	}
	
	if (accuracy <= 0.)
	{
		throw std::invalid_argument("Accuracy parameters must be greater than 0");
	}
	
	m_theta_th = theta;
	m_accuracy_th = accuracy;
}

void SimManager::launchSimulation()
{
	IncompressibleKernel flow_kernel(m_data, m_tMax, m_nb_steps, m_theta_inc, m_Lx, m_Ly, m_nx, m_ny, m_rho, m_mu, m_accuracy_u, m_accuracy_v, m_accuracy_p, m_clean_pressure);
	flow_kernel.defineBoundaryConditions(m_U_BC, m_V_BC, m_P_BC);
	flow_kernel.defineBody(m_body);
	flow_kernel.simulate_theta();

	ThermalKernel th_kernel(m_data, m_tMax, m_nb_steps, m_theta_th, m_accuracy_th, m_Lx, m_Ly, m_nx, m_ny, m_lambda, m_rho, m_cp, m_T_BC);
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