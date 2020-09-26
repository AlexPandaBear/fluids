#include "IncompressibleKernel.hxx"

IncompressibleKernel::IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu) :
	m_data(data),
	m_tMax(tMax),
	m_nb_steps(nb_steps),
	m_dt(tMax/nb_steps),
	m_Lx(Lx),
	m_Ly(Ly),
	m_nx(nx),
	m_ny(ny),
	m_dx(Lx/(nx-1)),
	m_dy(Ly/(ny-1)),
	m_rho(rho),
	m_mu(mu),
	m_nu(mu/rho)
{
	m_U = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));
	m_V = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));
	m_P = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));
}

IncompressibleKernel::~IncompressibleKernel() {}

void IncompressibleKernel::defineBoundaryConditions(std::vector<double> v_U_bc, std::vector<double> v_V_bc, std::vector<double> v_P_bc)
{
	size_t h;

	for (size_t i = 0; i < m_nx; i++)
	{
		m_U[i+1][0] = v_U_bc[i];
		m_V[i+1][0] = v_V_bc[i];
		m_P[i+1][0] = v_P_bc[i];
	}
	
	for (int i = 0; i < m_nx; i++)
	{
		m_U[i+1][m_ny+1] = v_U_bc[m_nx+i];
		m_V[i+1][m_ny+1] = v_V_bc[m_nx+i];
		m_P[i+1][m_ny+1] = v_P_bc[m_nx+i];
	}
	
	for (int i = 0; i < m_ny; i++)
	{
		m_U[0][i+1] = v_U_bc[2*m_nx+i];
		m_V[0][i+1] = v_V_bc[2*m_nx+i];
		m_P[0][i+1] = v_P_bc[2*m_nx+i];
	}
	
	for (int i = 0; i < m_ny; i++)
	{
		m_U[m_nx+1][i+1] = v_U_bc[2*m_nx+m_ny+i];
		m_V[m_nx+1][i+1] = v_V_bc[2*m_nx+m_ny+i];
		m_P[m_nx+1][i+1] = v_P_bc[2*m_nx+m_ny+i];
	}

	m_U[0][0] = 0.5 * (m_U[0][1] + m_U[1][0]);
	m_U[0][m_ny+1] = 0.5 * (m_U[0][m_ny] + m_U[1][m_ny+1]);
	m_U[m_nx+1][0] = 0.5 * (m_U[m_nx][0] + m_U[m_nx+1][1]);
	m_U[m_nx+1][m_ny+1] = 0.5 * (m_U[m_nx][m_ny+1] + m_U[m_nx+1][m_ny]);

	m_V[0][0] = 0.5 * (m_V[0][1] + m_V[1][0]);
	m_V[0][m_ny+1] = 0.5 * (m_V[0][m_ny] + m_V[1][m_ny+1]);
	m_V[m_nx+1][0] = 0.5 * (m_V[m_nx][0] + m_V[m_nx+1][1]);
	m_V[m_nx+1][m_ny+1] = 0.5 * (m_V[m_nx][m_ny+1] + m_V[m_nx+1][m_ny]);

	m_P[0][0] = 0.5 * (m_P[0][1] + m_P[1][0]);
	m_P[0][m_ny+1] = 0.5 * (m_P[0][m_ny] + m_P[1][m_ny+1]);
	m_P[m_nx+1][0] = 0.5 * (m_P[m_nx][0] + m_P[m_nx+1][1]);
	m_P[m_nx+1][m_ny+1] = 0.5 * (m_P[m_nx][m_ny+1] + m_P[m_nx+1][m_ny]);
}

void IncompressibleKernel::getPressureFieldAt(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_P[i+1][j+1] = m_data.getPressureAt(t, i, j);
		}
	}
}

void IncompressibleKernel::getVelocityFieldsAt(size_t t)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_U[i+1][j+1] = m_data.getXVelocityAt(t, i, j);
			m_V[i+1][j+1] = m_data.getYVelocityAt(t, i, j);
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


double IncompressibleKernel::d2udx2(size_t i, size_t j) const
{
	return (m_U[i+1][j] + m_U[i-1][j] - 2*m_U[i][j])/(m_dx*m_dx);
}

double IncompressibleKernel::d2udy2(size_t i, size_t j) const
{
	return (m_U[i][j+1] + m_U[i][j-1] - 2*m_U[i][j])/(m_dy*m_dy);
}

double IncompressibleKernel::d2vdx2(size_t i, size_t j) const
{
	return (m_V[i+1][j] + m_V[i-1][j] - 2*m_V[i][j])/(m_dx*m_dx);
}

double IncompressibleKernel::d2vdy2(size_t i, size_t j) const
{
	return (m_V[i][j+1] + m_V[i][j-1] - 2*m_V[i][j])/(m_dy*m_dy);
}


double IncompressibleKernel::d2uvdxdy(size_t i, size_t j) const
{
	return ((m_U[i+1][j+1]*m_V[i+1][j+1] - m_U[i+1][j-1]*m_V[i+1][j-1]) - (m_U[i-1][j+1]*m_V[i-1][j+1] - m_U[i-1][j-1]*m_V[i-1][j-1]))/(4*m_dx*m_dy);
}

double IncompressibleKernel::d2u2dx2(size_t i, size_t j) const
{
	return (pow(m_U[i+1][j], 2) + pow(m_U[i-1][j], 2) - 2*pow(m_U[i][j], 2))/(m_dx*m_dx);
}

double IncompressibleKernel::d2v2dy2(size_t i, size_t j) const
{
	return (pow(m_V[i][j+1], 2) + pow(m_V[i][j-1], 2) - 2*pow(m_V[i][j], 2))/(m_dy*m_dy);
}


double IncompressibleKernel::laplacian_u(size_t i, size_t j) const
{
	return d2udx2(i,j) + d2udy2(i,j);
}

double IncompressibleKernel::laplacian_v(size_t i, size_t j) const
{
	return d2vdx2(i,j) + d2vdy2(i,j);
}


double IncompressibleKernel::advection_u(size_t i, size_t j) const
{
	return m_U[i+1][j+1] * dudx(i, j) + m_V[i+1][j+1] * dudy(i, j);
}

double IncompressibleKernel::advection_v(size_t i, size_t j) const
{
	return m_U[i+1][j+1] * dvdx(i, j) + m_V[i+1][j+1] * dvdy(i, j);
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
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			m_data.setXVelocityAt(t+1, i, j, m_U[i+1][j+1] + m_dt * (diffusion_u(i+1, j+1) - advection_u(i+1, j+1) - pressure_u(i+1, j+1)));
			m_data.setYVelocityAt(t+1, i, j, m_V[i+1][j+1] + m_dt * (diffusion_v(i+1, j+1) - advection_v(i+1, j+1) - pressure_v(i+1, j+1)));
		}
	}
}

void IncompressibleKernel::computePressureField(PoissonSolver& solver, size_t t)
{
	std::vector<double> v_f(m_nx*m_ny, 0.);

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			v_f[i + m_nx*j] = - m_rho * (d2u2dx2(i+1,j+1) + 2.*d2uvdxdy(i+1,j+1) + d2v2dy2(i+1,j+1));
		}
	}

	solver.definePoissonFunction(v_f);
	solver.defineBoundaryConditions(0.);
	std::vector<double> v_P(solver.solve());

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			m_data.setPressureAt(t, i, j, v_P[i + m_nx*j]);
		}
	}
}

void IncompressibleKernel::printSimulationProgression(size_t t) const
{
	std::cout << "\rComputing step " << t+1 << " out of " << m_nb_steps << " -- " << 100.*(t+1)/m_nb_steps << "% completed     " << std::flush;
}

void IncompressibleKernel::simulate()
{
	PoissonSolver solver(m_Lx, m_Ly, m_nx, m_ny);
	std::cout << "Solving incompressible flow..." << std::endl;
	
	for (size_t t = 0; t < m_nb_steps; t++)
	{
		printSimulationProgression(t);
		getVelocityFieldsAt(t);
		computePressureField(solver, t);
		getPressureFieldAt(t);
		updateVelocityField(t);
	}

	std::cout << std::endl;
}