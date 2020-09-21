#include "IncompressibleKernel.hxx"

IncompressibleKernel::IncompressibleKernel() {}

IncompressibleKernel::~IncompressibleKernel() {}

void IncompressibleKernel::defineTimeParameters(double tMax, size_t nb_steps)
{
	m_tMax = tMax;
	m_nb_steps = nb_steps;
	m_dt = tMax/nb_steps;
}

void IncompressibleKernel::defineGridParameters(double Lx, double Ly, size_t nx, size_t ny)
{
	m_Lx = Lx;
	m_Ly = Ly;
	m_nx = nx;
	m_ny = ny;
	m_dx = Lx/(nx-1);
	m_dy = Ly/(ny-1);
}

void IncompressibleKernel::defineBodyShape(std::vector<std::vector<bool>> body_shape)
{
	v_body = body_shape;
}

void IncompressibleKernel::defineFluidProperties(double rho, double mu)
{
	m_rho = rho;
	m_mu = mu;
	m_nu = mu/rho;
}

void IncompressibleKernel::defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> P0)
{
	m_U = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	m_V = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	m_P = std::vector<std::vector<std::vector<double>>>(m_nb_steps+1, std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.)));
	m_U[0] = U0;
	m_V[0] = V0;
	m_P[0] = P0;
}

void IncompressibleKernel::initialize()
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			for (size_t t = 0; t < m_nb_steps+1; t++)
			{
				if (v_body[i][j])
				{
					m_U[t][i][j] = 0.;
					m_V[t][i][j] = 0.;
				}
			}
		}
	}
}

double IncompressibleKernel::dudx(size_t t, size_t i, size_t j) const
{
	return (m_U[t][i+1][j] - m_U[t][i-1][j])/(2*m_dx);
}

double IncompressibleKernel::dudy(size_t t, size_t i, size_t j) const
{
	return (m_U[t][i][j+1] - m_U[t][i][j-1])/(2*m_dy);
}

double IncompressibleKernel::dvdx(size_t t, size_t i, size_t j) const
{
	return (m_V[t][i+1][j] - m_V[t][i-1][j])/(2*m_dx);
}

double IncompressibleKernel::dvdy(size_t t, size_t i, size_t j) const
{
	return (m_V[t][i][j+1] - m_V[t][i][j-1])/(2*m_dy);
}

double IncompressibleKernel::laplacian_u(size_t t, size_t i, size_t j) const
{
	return (m_U[t][i+1][j] + m_U[t][i-1][j] - 2*m_U[t][i][j])/(m_dx*m_dx) + (m_U[t][i][j+1] + m_U[t][i][j-1] - 2*m_U[t][i][j])/(m_dy*m_dy);
}

double IncompressibleKernel::laplacian_v(size_t t, size_t i, size_t j) const
{
	return (m_V[t][i+1][j] + m_V[t][i-1][j] - 2*m_V[t][i][j])/(m_dx*m_dx) + (m_V[t][i][j+1] + m_V[t][i][j-1] - 2*m_V[t][i][j])/(m_dy*m_dy);
}

double IncompressibleKernel::advection_u(size_t t, size_t i, size_t j) const
{
	return m_U[t][i][j] * dudx(t, i, j) + m_V[t][i][j] * dudy(t, i, j);
}

double IncompressibleKernel::advection_v(size_t t, size_t i, size_t j) const
{
	return m_U[t][i][j] * dvdx(t, i, j) + m_V[t][i][j] * dvdy(t, i, j);
}

double IncompressibleKernel::diffusion_u(size_t t, size_t i, size_t j) const
{
	return m_nu * laplacian_u(t, i, j);
}

double IncompressibleKernel::diffusion_v(size_t t, size_t i, size_t j) const
{
	return m_nu * laplacian_v(t, i, j);
}

double IncompressibleKernel::pressure_u(size_t t, size_t i, size_t j) const
{
	return (m_P[t][i+1][j] - m_P[t][i-1][j])/(2*m_dx*m_rho);
}

double IncompressibleKernel::pressure_v(size_t t, size_t i, size_t j) const
{
	return (m_P[t][i][j+1] - m_P[t][i][j-1])/(2*m_dy*m_rho);
}

void IncompressibleKernel::updateVelocityField(size_t t)
{
	for (size_t i = 1; i < m_nx-1; i++)
	{
		for (size_t j = 1; j < m_ny-1; j++)
		{
			m_U[t+1][i][j] = m_U[t][i][j] + m_dt * (diffusion_u(t, i, j) - advection_u(t, i, j) - pressure_u(t, i, j));
			m_V[t+1][i][j] = m_V[t][i][j] + m_dt * (diffusion_v(t, i, j) - advection_v(t, i, j) - pressure_v(t, i, j));
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		m_U[t+1][0][j] = m_U[t][0][j];
		m_U[t+1][m_nx-1][j] = m_U[t][m_nx-1][j];
		m_V[t+1][0][j] = m_V[t][0][j];
		m_V[t+1][m_nx-1][j] = m_V[t][m_nx-1][j];
	}

	for (size_t i = 0; i < m_nx; i++)
	{
		m_U[t+1][i][0] = m_U[t][i][0];
		m_U[t+1][i][m_ny-1] = m_U[t][i][m_ny-1];
		m_V[t+1][i][0] = m_V[t][i][0];
		m_V[t+1][i][m_ny-1] = m_V[t][i][m_ny-1];
	}
}

void IncompressibleKernel::updatePressureField(size_t t, double mean_error_max)
{
	double mean_err = mean_error_max+1;

	double coef_1 = 0.5*m_dy*m_dy/(m_dx*m_dx + m_dy*m_dy);
	double coef_2 = 0.5*m_dx*m_dx/(m_dx*m_dx + m_dy*m_dy);
	double coef_3 = 0.5*m_rho*m_dx*m_dx*m_dy*m_dy/(m_dx*m_dx + m_dy*m_dy);

	std::vector<std::vector<double>> v_P_tmp = m_P[t];

	while (mean_err > mean_error_max)
	{
		mean_err = 0.;

		for (size_t i = 1; i < m_nx-1; i++)
		{
			for (size_t j = 1; j < m_ny-1; j++)
			{
				m_P[t+1][i][j] = coef_1*(v_P_tmp[i+1][j] + v_P_tmp[i-1][j]) + coef_2*(v_P_tmp[i][j+1] + v_P_tmp[i][j-1]) + coef_3*(pow(dudx(t, i, j), 2) + 2.*dudy(t, i, j)*dvdx(t, i, j) + pow(dvdy(t, i, j), 2));
				mean_err += abs(m_P[t+1][i][j] - v_P_tmp[i][j]);
				v_P_tmp[i][j] = m_P[t+1][i][j];
			}
		}

		mean_err /= m_nx*m_ny;

		for (size_t j = 0; j < m_ny; j++)
		{
			m_P[t+1][0][j] = m_P[t][0][j];
			m_P[t+1][m_nx-1][j] = m_P[t][m_nx-1][j];
		}

		for (size_t i = 0; i < m_nx; i++)
		{
			m_P[t+1][i][0] = m_P[t][i][0];
			m_P[t+1][i][m_ny-1] = m_P[t][i][m_ny-1];
		}
	}
}

void IncompressibleKernel::simulate()
{
	for (size_t t = 0; t < m_nb_steps; t++)
	{
		updateVelocityField(t);
		updatePressureField(t, 100.);
	}
}

std::vector<std::vector<std::vector<double>>> IncompressibleKernel::getU() const
{
	return m_U;
}

std::vector<std::vector<std::vector<double>>> IncompressibleKernel::getV() const
{
	return m_V;
}

std::vector<std::vector<std::vector<double>>> IncompressibleKernel::getP() const
{
	return m_P;
}