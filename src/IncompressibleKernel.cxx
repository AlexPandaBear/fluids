#include "IncompressibleKernel.hxx"

IncompressibleKernel::IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu, double err_max) :
	m_data(data),
	m_tMax(tMax),
	m_nb_steps(nb_steps),
	m_Lx(Lx),
	m_Ly(Ly),
	m_nx(nx),
	m_ny(ny),
	m_rho(rho),
	m_mu(mu),
	m_err_max(err_max)
{
	m_dt = m_tMax/m_nb_steps;
	m_dx = m_Lx/(m_nx-1);
	m_dy = m_Ly/(m_ny-1);

	m_nu = m_mu/m_rho;

	m_U = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_V = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_P = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_P_tmp = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));

	m_coef_1 = 0.5*m_dy*m_dy/(m_dx*m_dx + m_dy*m_dy);
	m_coef_2 = 0.5*m_dx*m_dx/(m_dx*m_dx + m_dy*m_dy);
	m_coef_3 = 0.5*m_rho*m_dx*m_dx*m_dy*m_dy/(m_dx*m_dx + m_dy*m_dy);
}

IncompressibleKernel::~IncompressibleKernel() {}

void IncompressibleKernel::initialize()
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			if (v_body[i][j])
			{
				m_U[i][j] = 0.;
				m_V[i][j] = 0.;
			}
		}
	}
}

void IncompressibleKernel::getAllFieldsAt(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_U[i][j] = m_data.getXVelocityAt(t, i, j);
			m_V[i][j] = m_data.getYVelocityAt(t, i, j);
			m_P[i][j] = m_data.getPressureAt(t, i, j);
		}
	}
}

void IncompressibleKernel::getVelocityFieldsAt(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_U[i][j] = m_data.getXVelocityAt(t, i, j);
			m_V[i][j] = m_data.getYVelocityAt(t, i, j);
		}
	}
}

double IncompressibleKernel::dudx(size_t i, size_t j) const
{
	return (m_U[i+1][j] - m_U[i-1][j])/(2*m_dx);
}

double IncompressibleKernel::dudy(size_t i, size_t j) const
{
	return (m_U[i][j+1] - m_U[i][j-1])/(2*m_dy);
}

double IncompressibleKernel::dvdx(size_t i, size_t j) const
{
	return (m_V[i+1][j] - m_V[i-1][j])/(2*m_dx);
}

double IncompressibleKernel::dvdy(size_t i, size_t j) const
{
	return (m_V[i][j+1] - m_V[i][j-1])/(2*m_dy);
}

double IncompressibleKernel::laplacian_u(size_t i, size_t j) const
{
	return (m_U[i+1][j] + m_U[i-1][j] - 2*m_U[i][j])/(m_dx*m_dx) + (m_U[i][j+1] + m_U[i][j-1] - 2*m_U[i][j])/(m_dy*m_dy);
}

double IncompressibleKernel::laplacian_v(size_t i, size_t j) const
{
	return (m_V[i+1][j] + m_V[i-1][j] - 2*m_V[i][j])/(m_dx*m_dx) + (m_V[i][j+1] + m_V[i][j-1] - 2*m_V[i][j])/(m_dy*m_dy);
}

double IncompressibleKernel::advection_u(size_t i, size_t j) const
{
	return m_U[i][j] * dudx(i, j) + m_V[i][j] * dudy(i, j);
}

double IncompressibleKernel::advection_v(size_t i, size_t j) const
{
	return m_U[i][j] * dvdx(i, j) + m_V[i][j] * dvdy(i, j);
}

double IncompressibleKernel::diffusion_u(size_t i, size_t j) const
{
	return m_nu * laplacian_u(i, j);
}

double IncompressibleKernel::diffusion_v(size_t i, size_t j) const
{
	return m_nu * laplacian_v(i, j);
}

double IncompressibleKernel::pressure_u(size_t i, size_t j) const
{
	return (m_P[i+1][j] - m_P[i-1][j])/(2*m_dx*m_rho);
}

double IncompressibleKernel::pressure_v(size_t i, size_t j) const
{
	return (m_P[i][j+1] - m_P[i][j-1])/(2*m_dy*m_rho);
}

void IncompressibleKernel::updateVelocityField(size_t t)
{
	for (size_t i = 1; i < m_nx-1; i++)
	{
		for (size_t j = 1; j < m_ny-1; j++)
		{
			m_data.setXVelocityAt(t+1, i, j, m_U[i][j] + m_dt * (diffusion_u(i, j) - advection_u(i, j) - pressure_u(i, j)));
			m_data.setYVelocityAt(t+1, i, j, m_V[i][j] + m_dt * (diffusion_v(i, j) - advection_v(i, j) - pressure_v(i, j)));
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		m_data.setXVelocityAt(t+1, 0, j, m_U[0][j]);
		m_data.setXVelocityAt(t+1, m_nx-1, j, m_U[m_nx-1][j]);
		m_data.setYVelocityAt(t+1, 0, j, m_V[0][j]);
		m_data.setYVelocityAt(t+1, m_nx-1, j, m_V[m_nx-1][j]);
	}

	for (size_t i = 0; i < m_nx; i++)
	{
		m_data.setXVelocityAt(t+1, i, 0, m_U[i][0]);
		m_data.setXVelocityAt(t+1, i, m_ny-1, m_U[i][m_ny-1]);
		m_data.setYVelocityAt(t+1, i, 0, m_V[i][0]);
		m_data.setYVelocityAt(t+1, i, m_ny-1, m_V[i][m_ny-1]);
	}
}

void IncompressibleKernel::updatePressureField(size_t t, double mean_error_max)
{
	double mean_err = mean_error_max+1;

	while (mean_err > mean_error_max)
	{
		for (size_t i = 1; i < m_nx-1; i++)
		{
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_P_tmp[i][j] = m_coef_1*(m_P[i+1][j] + m_P[i-1][j]) + m_coef_2*(m_P[i][j+1] + m_P[i][j-1]) + m_coef_3*(pow(dudx(i, j), 2) + 2.*dudy(i, j)*dvdx(i, j) + pow(dvdy(i, j), 2));
			}
		}

		for (size_t j = 0; j < m_ny; j++)
		{
			m_P_tmp[0][j] = m_P[0][j];
			m_P_tmp[m_nx-1][j] = m_P[m_nx-1][j];
		}

		for (size_t i = 0; i < m_nx; i++)
		{
			m_P_tmp[i][0] = m_P[i][0];
			m_P_tmp[i][m_ny-1] = m_P[i][m_ny-1];
		}

		mean_err = 0.;

		for (size_t i = 0; i < m_nx; i++)
		{
			for (size_t j = 0; j < m_ny; j++)
			{
				mean_err += m_P_tmp[i][j] - m_P[i][j];
				m_P[i][j] = m_P_tmp[i][j];
			}
		}

		mean_err /= m_nx*m_ny;
	}

	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_data.setPressureAt(t+1, i, j, m_P[i][j]);
		}
	}
}

void IncompressibleKernel::printSimulationProgression(size_t t) const
{
	std::cout << "\rComputing step " << t+1 << " out of " << m_nb_steps << " -- " << 100.*(t+1)/m_nb_steps << "% completed     " << std::flush;
}

void IncompressibleKernel::simulate()
{
	std::cout << "Solving incompressible flow..." << std::endl;

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		printSimulationProgression(t);
		getAllFieldsAt(t);
		updateVelocityField(t);
		getVelocityFieldsAt(t+1);
		updatePressureField(t, m_err_max);
	}

	std::cout << std::endl;
}