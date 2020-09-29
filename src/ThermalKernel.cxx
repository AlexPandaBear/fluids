#include "ThermalKernel.hxx"

ThermalKernel::ThermalKernel(DataKeeper& data, double tMax, size_t nb_steps, double theta, double accuracy, double Lx, double Ly, size_t nx, size_t ny, double lambda, double rho, double cp, double T_BC) :
	m_data(data),
	m_tMax(tMax),
	m_nb_steps(nb_steps),
	m_dt(tMax/nb_steps),
	m_theta(theta),
	m_accuracy(accuracy),
	m_Lx(Lx),
	m_Ly(Ly),
	m_nx(nx),
	m_ny(ny),
	m_dx(Lx/(nx-1)),
	m_dy(Ly/(ny-1)),
	m_nb_eq(nx*ny),
	m_lambda(lambda),
	m_rho(rho),
	m_cp(cp),
	m_K(lambda/(rho*cp)),
	m_T_BC(T_BC),
	m_kx(m_K*m_dt/(m_dx*m_dx)),
	m_ky(m_K*m_dt/(m_dy*m_dy))
{
	m_u = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_v = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	m_T = std::vector<std::vector<double>>(m_nx, std::vector<double>(m_ny, 0.));
	
	mat_A_a = std::vector<double>(m_nb_eq, 1. + 2.*m_theta*(m_kx + m_ky));
	mat_A_b = std::vector<double>(m_nb_eq, 0.);
	mat_A_c = std::vector<double>(m_nb_eq, 0.);
	mat_A_d = std::vector<double>(m_nb_eq, 0.);
	mat_A_e = std::vector<double>(m_nb_eq, 0.);

	vec_X = std::vector<double>(m_nb_eq, 0.);
	vec_B = std::vector<double>(m_nb_eq, 0.);
}

ThermalKernel::~ThermalKernel() {}

void ThermalKernel::display()
{
	for (size_t l = 0; l < m_nb_eq; l++)
	{
		std::cout << mat_A_e[l] << "..." << mat_A_c[l] << " " << mat_A_a[l] << " " << mat_A_b[l] << "..." << mat_A_d[l] << std::endl;
	}
}

void ThermalKernel::computeMatrixA(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_u[i][j] = (0.5*m_dt/m_dx)*m_data.getXVelocityAt(t+1, i, j);
			m_v[i][j] = (0.5*m_dt/m_dy)*m_data.getYVelocityAt(t+1, i, j);
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			mat_A_b[i + m_nx*j] = - m_theta * (m_kx - m_u[i][j]);
			mat_A_c[i + m_nx*j] = - m_theta * (m_kx + m_u[i][j]);
			mat_A_d[i + m_nx*j] = - m_theta * (m_ky - m_v[i][j]);
			mat_A_e[i + m_nx*j] = - m_theta * (m_ky + m_v[i][j]);
		}
	}
}

void ThermalKernel::computeVectorB(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_u[i][j] = (0.5*m_dt/m_dx)*m_data.getXVelocityAt(t, i, j);
			m_v[i][j] = (0.5*m_dt/m_dy)*m_data.getYVelocityAt(t, i, j);
			m_T[i][j] = m_data.getTemperatureAt(t, i, j);

			vec_B[i + m_nx*j] = (1. - 2.*(1. - m_theta)*(m_kx + m_ky)) * m_T[i][j];
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		vec_B[m_nx*j] += (1. - m_theta)*(m_kx - m_u[0][j]) * m_T[1][j];
		vec_B[m_nx*j] += (1. - m_theta)*(m_kx + m_u[0][j]) * m_T_BC;

		for (size_t i = 1; i < m_nx-1; i++)
		{
			vec_B[i + m_nx*j] += (1. - m_theta)*(m_kx - m_u[i][j]) * m_T[i+1][j];
			vec_B[i + m_nx*j] += (1. - m_theta)*(m_kx + m_u[i][j]) * m_T[i-1][j];
		}

		vec_B[(m_nx-1) + m_nx*j] += (1. - m_theta)*(m_kx - m_u[m_nx-1][j]) * m_T_BC;
		vec_B[(m_nx-1) + m_nx*j] += (1. - m_theta)*(m_kx + m_u[m_nx-1][j]) * m_T[m_nx-2][j];
	}

	for (size_t i = 0; i < m_nx; i++)
	{
		vec_B[i] += (1. - m_theta)*(m_ky - m_v[i][0]) * m_T[i][1];
		vec_B[i] += (1. - m_theta)*(m_ky + m_v[i][0]) * m_T_BC;

		for (size_t j = 1; j < m_ny-1; j++)
		{
			vec_B[i + m_nx*j] += (1. - m_theta)*(m_ky - m_v[i][j]) * m_T[i][j+1];
			vec_B[i + m_nx*j] += (1. - m_theta)*(m_ky + m_v[i][j]) * m_T[i][j-1];
		}

		vec_B[i + m_nx*(m_ny-1)] += (1. - m_theta)*(m_ky - m_v[i][m_ny-1]) * m_T_BC;
		vec_B[i + m_nx*(m_ny-1)] += (1. - m_theta)*(m_ky + m_v[i][m_ny-1]) * m_T[i][m_ny-2];
	}
}

void ThermalKernel::solve_GaussSeidel()
{
	double mean_diff;
	double new_xi;

	do
	{
		mean_diff = 0.;

		for (size_t i = 0; i < m_nx; i++)
		{
			for (size_t j = 0; j < m_ny; j++)
			{
				new_xi = vec_B[i + m_nx*j];

				if (i == m_nx-1)
				{
					new_xi -= mat_A_b[i + m_nx*j] * m_T_BC;
					//new_xi -= mat_A_b[i + m_nx*j] * vec_X[m_ny*j];
				}
				else
				{
					new_xi -= mat_A_b[i + m_nx*j] * vec_X[i + m_nx*j + 1];
				}

				if (i == 0)
				{
					new_xi -= mat_A_c[i + m_nx*j] * m_T_BC;
					//new_xi -= mat_A_c[i + m_nx*j] * vec_X[m_nx*(j+1)-1];
				}
				else
				{
					new_xi -= mat_A_c[i + m_nx*j] * vec_X[i + m_nx*j - 1];
				}

				if (j == m_ny-1)
				{
					new_xi -= mat_A_d[i + m_nx*j] * m_T_BC;
				}
				else
				{
					new_xi -= mat_A_d[i + m_nx*j] * vec_X[i + m_nx*j + m_nx];
				}

				if (j == 0)
				{
					new_xi -= mat_A_e[i + m_nx*j] * m_T_BC;
				}
				else
				{
					new_xi -= mat_A_e[i + m_nx*j] * vec_X[i + m_nx*j - m_nx];
				}
				
				new_xi /= mat_A_a[i + m_nx*j];
				mean_diff += abs(new_xi - vec_X[i + m_nx*j]);
				vec_X[i + m_nx*j] = new_xi;
			}

		}

		mean_diff /= m_nb_eq;
		
	}
	while (mean_diff > m_accuracy);
}

void ThermalKernel::printSimulationProgression(size_t t) const
{
	std::cout << "\rComputing step " << t+1 << " out of " << m_nb_steps << " -- " << 100.*(t+1)/m_nb_steps << "% completed     " << std::flush;
}
	
void ThermalKernel::simulate()
{
	std::cout << "Solving heat equation..." << std::endl;

	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			vec_X[i + m_nx*j] = m_data.getTemperatureAt(0, i, j);
		}
	}

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		printSimulationProgression(t);
		computeMatrixA(t);
		computeVectorB(t);
		//display();
		solve_GaussSeidel();

		for (size_t j = 0; j < m_ny; j++)
		{
			//std::cout << std::endl;
			for (size_t i = 0; i < m_nx; i++)
			{
				m_data.setTemperatureAt(t+1, i, j, vec_X[i + m_nx*j]);
				//std::cout << vec_X[i + m_nx*j] << "  ";
			}
		}
	}

	std::cout << std::endl;
}