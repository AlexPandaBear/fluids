#include "ThermalKernel.hxx"

ThermalKernel::ThermalKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double lambda, double rho, double cv, bool temp_BC, bool flux_BC, double BC_value) :
	m_data(data),
	m_tMax(tMax),
	m_nb_steps(nb_steps),
	m_Lx(Lx),
	m_Ly(Ly),
	m_nx(nx),
	m_ny(ny),
	m_lambda(lambda),
	m_rho(rho),
	m_cv(cv),
	m_temp_BC(temp_BC),
	m_flux_BC(flux_BC),
	m_BC_value(BC_value)
{
	m_dt = m_tMax/m_nb_steps;
	m_dx = m_Lx/(m_nx-1);
	m_dy = m_Ly/(m_ny-1);

	m_Cv = m_rho*m_cv*m_dx*m_dy;

	m_coef_x = m_lambda*m_dy/m_Cv;
	m_coef_y = m_lambda*m_dx/m_Cv;

	m_T = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_U = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_V = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
}

ThermalKernel::~ThermalKernel() {}

void ThermalKernel::getAllFieldsAt(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_T[i][j] = m_data.getTemperatureAt(t, i, j);
			m_U[i][j] = m_data.getXVelocityAt(t, i, j);
			m_V[i][j] = m_data.getYVelocityAt(t, i, j);
		}
	}
}

double ThermalKernel::diffusion(size_t i, size_t j) const
{
	return m_coef_x*(m_T[i+1][j] + m_T[i-1][j] - 2*m_T[i][j]) + m_coef_y*(m_T[i][j+1] + m_T[i][j-1] - 2*m_T[i][j]);
}

double ThermalKernel::advection(size_t i, size_t j) const
{
	return m_U[i][j]*(m_T[i+1][j]-m_T[i-1][j])/m_dx + m_V[i][j]*(m_T[i][j+1]-m_T[i][j-1])/m_dy;
}

double ThermalKernel::creation(size_t i, size_t j) const
{
	//return m_nu*((m_U[i+1][j] + m_U[i-1][j] - 2*m_U[i][j])/(m_dy*m_dy) + (m_V[i+1][j] + m_V[i-1][j] - 2*m_V[i][j])/(m_dx*m_dx));
	return 0.;
}
	
void ThermalKernel::simulate()
{
	for (size_t t = 0; t < m_nb_steps; t++)
	{
		getAllFieldsAt(t);

		for (size_t i = 1; i < m_nx-1; i++)
		{
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_data.setTemperatureAt(t+1, i, j, m_T[i][j] + m_dt * (diffusion(i, j) - advection(i, j) + creation(i, j)));
			}
		}

		if (m_temp_BC)
		{
			for (size_t i = 0; i < m_nx; i++)
			{
				m_data.setTemperatureAt(t+1, i, 0, m_BC_value);
				m_data.setTemperatureAt(t+1, i, m_ny-1, m_BC_value);
			}
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_data.setTemperatureAt(t+1, 0, j, m_BC_value);
				m_data.setTemperatureAt(t+1, m_nx-1, j, m_BC_value);
			}
		}
		else if (m_flux_BC)
		{
			for (size_t i = 1; i < m_nx-1; i++)
			{
				m_data.setTemperatureAt(t+1, i, 0, m_T[i][0] + m_dt * (m_coef_x*(m_T[i+1][0] + m_T[i-1][0] - 2*m_T[i][0]) + m_coef_y*m_BC_value));
				m_data.setTemperatureAt(t+1, i, m_ny-1, m_T[i][m_ny-1] + m_dt * (m_coef_x*(m_T[i+1][m_ny-1] + m_T[i-1][m_ny-1] - 2*m_T[i][m_ny-1]) + m_coef_y*m_BC_value));
			}
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_data.setTemperatureAt(t+1, 0, j, m_T[0][j] + m_dt * (m_coef_y*(m_T[0][j+1] + m_T[0][j-1] - 2*m_T[0][j]) + m_coef_x*m_BC_value));
				m_data.setTemperatureAt(t+1, m_nx-1, j, m_T[m_nx-1][j] + m_dt * (m_coef_y*(m_T[m_nx-1][j+1] + m_T[m_nx-1][j-1] - 2*m_T[m_nx-1][j]) + m_coef_x*m_BC_value));
			}
			m_data.setTemperatureAt(t+1, 0, 0, m_T[0][0] + m_dt * (m_coef_x + m_coef_y) * m_BC_value);
			m_data.setTemperatureAt(t+1, 0, m_ny-1, m_T[0][m_ny-1] + m_dt * (m_coef_x + m_coef_y) * m_BC_value);
			m_data.setTemperatureAt(t+1, m_nx-1, 0, m_T[m_nx-1][0] + m_dt * (m_coef_x + m_coef_y) * m_BC_value);
			m_data.setTemperatureAt(t+1, m_nx-1, m_ny-1, m_T[m_nx-1][m_ny-1] + m_dt * (m_coef_x + m_coef_y) * m_BC_value);
		}
		else
		{
			throw std::logic_error("Bad definition of boundary conditions");
		}
	}
}