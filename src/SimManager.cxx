#include "SimManager.hxx"

SimManager::SimManager() {}
SimManager::~SimManager() {}

double SimManager::diffusion(size_t t, size_t i, size_t j) const
{
	return m_coef_x*(m_T[t][i+1][j] + m_T[t][i-1][j] - 2*m_T[t][i][j]) + m_coef_y*(m_T[t][i][j+1] + m_T[t][i][j-1] - 2*m_T[t][i][j]);
}

double SimManager::advection(size_t t, size_t i, size_t j) const
{
	return m_U[t][i][j]*(m_T[t][i+1][j]-m_T[t][i-1][j])/m_dx + m_V[t][i][j]*(m_T[t][i][j+1]-m_T[t][i][j-1])/m_dy;
}

double SimManager::creation(size_t t, size_t i, size_t j) const
{
	//return m_nu*((m_U[t][i+1][j] + m_U[t][i-1][j] - 2*m_U[t][i][j])/(m_dy*m_dy) + (m_V[t][i+1][j] + m_V[t][i-1][j] - 2*m_V[t][i][j])/(m_dx*m_dx));
	return 0.;
}

void SimManager::defineTimeParameters(double tMax, size_t nb_steps)
{
	m_tMax = tMax;
	m_nb_steps = nb_steps;
	m_dt = tMax/nb_steps;
}

void SimManager::defineGridParameters(double Lx, double Ly, size_t nx, size_t ny)
{
	m_Lx = Lx;
	m_Ly = Ly;
	m_nx = nx;
	m_ny = ny;
	m_dx = Lx/(nx-1);
	m_dy = Ly/(ny-1);
}

void SimManager::defineFluidProperties(double lambda, double rho, double cv, double mu)
{
	m_lambda = lambda;
	m_rho = rho;
	m_cv = cv;
	m_mu = mu;
	m_nu = mu/rho;
	m_Cv = rho*cv*m_dx*m_dy;
	m_coef_x = lambda*m_dy/m_Cv;
	m_coef_y = lambda*m_dx/m_Cv;
}

void SimManager::defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> P0, std::vector<std::vector<double>> T0)
{
	m_U = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	m_V = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	m_P = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	m_T = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	
	m_U[0] = U0;
	m_V[0] = V0;
	m_P[0] = P0;
	m_T[0] = T0;
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

void SimManager::computeFlow()
{
	IncompressibleKernel kernel;
	
	kernel.defineTimeParameters(m_tMax, m_nb_steps);
	kernel.defineGridParameters(m_Lx, m_Ly, m_nx, m_ny);
	//kernel.defineBodyShape(...)
	kernel.defineFluidProperties(m_rho, m_mu);
	kernel.defineInitialState(m_U[0], m_V[0], m_P[0]);
	
	kernel.simulate();

	m_U = kernel.getU();
	m_V = kernel.getV();
	m_P = kernel.getP();
}

void SimManager::launchSimulation()
{
	computeFlow();

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		for (size_t i = 1; i < m_nx-1; i++)
		{
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_T[t+1][i][j] = m_T[t][i][j] + m_dt * (diffusion(t, i, j) - advection(t, i, j) + creation(t, i, j));
			}
		}

		if (m_temp_BC)
		{
			for (size_t i = 0; i < m_nx; i++)
			{
				m_T[t+1][i][0] = m_BC_value;
				m_T[t+1][i][m_ny-1] = m_BC_value;
			}
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_T[t+1][0][j] = m_BC_value;
				m_T[t+1][m_nx-1][j] = m_BC_value;
			}
		}
		else if (m_flux_BC)
		{
			for (size_t i = 1; i < m_nx-1; i++)
			{
				m_T[t+1][i][0] = m_T[t][i][0] + m_dt * (m_coef_x*(m_T[t][i+1][0] + m_T[t][i-1][0] - 2*m_T[t][i][0]) + m_coef_y*m_BC_value);
				m_T[t+1][i][m_ny-1] = m_T[t][i][m_ny-1] + m_dt * (m_coef_x*(m_T[t][i+1][m_ny-1] + m_T[t][i-1][m_ny-1] - 2*m_T[t][i][m_ny-1]) + m_coef_y*m_BC_value);
			}
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_T[t+1][0][j] = m_T[t][0][j] + m_dt * (m_coef_y*(m_T[t][0][j+1] + m_T[t][0][j-1] - 2*m_T[t][0][j]) + m_coef_x*m_BC_value);
				m_T[t+1][m_nx-1][j] = m_T[t][m_nx-1][j] + m_dt * (m_coef_y*(m_T[t][m_nx-1][j+1] + m_T[t][m_nx-1][j-1] - 2*m_T[t][m_nx-1][j]) + m_coef_x*m_BC_value);
			}
			m_T[t+1][0][0] = m_T[t][0][0] + m_dt * (m_coef_x + m_coef_y) * m_BC_value;
			m_T[t+1][0][m_ny-1] = m_T[t][0][m_ny-1] + m_dt * (m_coef_x + m_coef_y) * m_BC_value;
			m_T[t+1][m_nx-1][0] = m_T[t][m_nx-1][0] + m_dt * (m_coef_x + m_coef_y) * m_BC_value;
			m_T[t+1][m_nx-1][m_ny-1] = m_T[t][m_nx-1][m_ny-1] + m_dt * (m_coef_x + m_coef_y) * m_BC_value;
		}
		else
		{
			throw std::logic_error("Bad definition of boundary conditions");
		}
	}
}

std::vector<std::vector<std::vector<double>>> SimManager::getResults() const
{
	return m_T;
}