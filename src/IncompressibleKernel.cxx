#include "IncompressibleKernel.hxx"

IncompressibleKernel::IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double theta, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu, double accuracy_u, double accuracy_v, double accuracy_p, bool clean_pressure) :
	m_data(data),
	m_tMax(tMax),
	m_nb_steps(nb_steps),
	m_dt(tMax/nb_steps),
	m_theta(theta),
	m_Lx(Lx),
	m_Ly(Ly),
	m_nx(nx),
	m_ny(ny),
	m_nb_pts(nx*ny),
	m_dx(Lx/(nx-1)),
	m_dy(Ly/(ny-1)),
	m_dxdy(m_dx/m_dy),
	m_dtdx(m_dt/m_dx),
	m_dtdy(m_dt/m_dy),
	m_rho(rho),
	m_mu(mu),
	m_nu(mu/rho),
	m_accuracy_u(accuracy_u),
	m_accuracy_v(accuracy_v),
	m_accuracy_p(accuracy_p),
	m_clean_pressure(clean_pressure),
	m_nu_x(m_nu*m_dt/pow(m_dx, 2)),
	m_nu_y(m_nu*m_dt/pow(m_dy, 2)),
	m_delta_x(0.5*m_dt/(m_rho*m_dx)),
	m_delta_y(0.5*m_dt/(m_rho*m_dy))
{
	m_U = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));
	m_V = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));
	m_P = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));

	U_tilde = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));
	V_tilde = std::vector<std::vector<double>>(m_nx+2, std::vector<double>(m_ny+2, 0.));

	diag_alpha = std::vector<double>(m_nb_pts, 1. + 2.*m_theta*(m_nu_x + m_nu_y));
	diag_beta1 = std::vector<double>(m_nb_pts, 0.);
	diag_beta2 = std::vector<double>(m_nb_pts, 0.);
	diag_gamma1 = std::vector<double>(m_nb_pts, 0.);
	diag_gamma2 = std::vector<double>(m_nb_pts, 0.);

	diag_delta1u = std::vector<double>(m_nb_pts, m_delta_x);
	diag_delta2u = std::vector<double>(m_nb_pts, -m_delta_x);
	diag_delta1v = std::vector<double>(m_nb_pts, m_delta_y);
	diag_delta2v = std::vector<double>(m_nb_pts, -m_delta_y);

	diag_Au = std::vector<double>(m_nb_pts, 0.);
	diag_B1u = std::vector<double>(m_nb_pts, 0.);
	diag_B2u = std::vector<double>(m_nb_pts, 0.);
	diag_D1 = std::vector<double>(m_nb_pts, 0.);
	diag_D2 = std::vector<double>(m_nb_pts, 0.);
	diag_D3 = std::vector<double>(m_nb_pts, 0.);
	diag_D4 = std::vector<double>(m_nb_pts, 0.);

	diag_Av = std::vector<double>(m_nb_pts, 0.);
	diag_C1v = std::vector<double>(m_nb_pts, 0.);
	diag_C2v = std::vector<double>(m_nb_pts, 0.);
	diag_E1 = std::vector<double>(m_nb_pts, 0.);
	diag_E2 = std::vector<double>(m_nb_pts, 0.);
	diag_E3 = std::vector<double>(m_nb_pts, 0.);
	diag_E4 = std::vector<double>(m_nb_pts, 0.);

	diag_A = std::vector<double>(m_nb_pts, -2.*(1. + pow(m_dxdy, 2)));
	diag_B1 = std::vector<double>(m_nb_pts, 1.);
	diag_B2 = std::vector<double>(m_nb_pts, 1.);
	diag_C1 = std::vector<double>(m_nb_pts, pow(m_dxdy, 2));
	diag_C2 = std::vector<double>(m_nb_pts, pow(m_dxdy, 2));

	vec_X = std::vector<double>(3*m_nb_pts, 0.);
	vec_B = std::vector<double>(3*m_nb_pts, 0.);
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

		U_tilde[i+1][0] = v_U_bc[i];
		V_tilde[i+1][0] = v_V_bc[i];
	}
	
	for (int i = 0; i < m_nx; i++)
	{
		m_U[i+1][m_ny+1] = v_U_bc[m_nx+i];
		m_V[i+1][m_ny+1] = v_V_bc[m_nx+i];
		m_P[i+1][m_ny+1] = v_P_bc[m_nx+i];

		U_tilde[i+1][m_ny+1] = v_U_bc[m_nx+i];
		V_tilde[i+1][m_ny+1] = v_V_bc[m_nx+i];
	}
	
	for (int i = 0; i < m_ny; i++)
	{
		m_U[0][i+1] = v_U_bc[2*m_nx+i];
		m_V[0][i+1] = v_V_bc[2*m_nx+i];
		m_P[0][i+1] = v_P_bc[2*m_nx+i];

		U_tilde[0][i+1] = v_U_bc[2*m_nx+i];
		V_tilde[0][i+1] = v_V_bc[2*m_nx+i];
	}
	
	for (int i = 0; i < m_ny; i++)
	{
		m_U[m_nx+1][i+1] = v_U_bc[2*m_nx+m_ny+i];
		m_V[m_nx+1][i+1] = v_V_bc[2*m_nx+m_ny+i];
		m_P[m_nx+1][i+1] = v_P_bc[2*m_nx+m_ny+i];

		U_tilde[m_nx+1][i+1] = v_U_bc[2*m_nx+m_ny+i];
		V_tilde[m_nx+1][i+1] = v_V_bc[2*m_nx+m_ny+i];
	}

	m_U[0][0] = 0.5 * (m_U[0][1] + m_U[1][0]);
	m_U[0][m_ny+1] = 0.5 * (m_U[0][m_ny] + m_U[1][m_ny+1]);
	m_U[m_nx+1][0] = 0.5 * (m_U[m_nx][0] + m_U[m_nx+1][1]);
	m_U[m_nx+1][m_ny+1] = 0.5 * (m_U[m_nx][m_ny+1] + m_U[m_nx+1][m_ny]);

	U_tilde[0][0] = 0.5 * (m_U[0][1] + m_U[1][0]);
	U_tilde[0][m_ny+1] = 0.5 * (m_U[0][m_ny] + m_U[1][m_ny+1]);
	U_tilde[m_nx+1][0] = 0.5 * (m_U[m_nx][0] + m_U[m_nx+1][1]);
	U_tilde[m_nx+1][m_ny+1] = 0.5 * (m_U[m_nx][m_ny+1] + m_U[m_nx+1][m_ny]);

	m_V[0][0] = 0.5 * (m_V[0][1] + m_V[1][0]);
	m_V[0][m_ny+1] = 0.5 * (m_V[0][m_ny] + m_V[1][m_ny+1]);
	m_V[m_nx+1][0] = 0.5 * (m_V[m_nx][0] + m_V[m_nx+1][1]);
	m_V[m_nx+1][m_ny+1] = 0.5 * (m_V[m_nx][m_ny+1] + m_V[m_nx+1][m_ny]);

	V_tilde[0][0] = 0.5 * (m_V[0][1] + m_V[1][0]);
	V_tilde[0][m_ny+1] = 0.5 * (m_V[0][m_ny] + m_V[1][m_ny+1]);
	V_tilde[m_nx+1][0] = 0.5 * (m_V[m_nx][0] + m_V[m_nx+1][1]);
	V_tilde[m_nx+1][m_ny+1] = 0.5 * (m_V[m_nx][m_ny+1] + m_V[m_nx+1][m_ny]);

	m_P[0][0] = 0.5 * (m_P[0][1] + m_P[1][0]);
	m_P[0][m_ny+1] = 0.5 * (m_P[0][m_ny] + m_P[1][m_ny+1]);
	m_P[m_nx+1][0] = 0.5 * (m_P[m_nx][0] + m_P[m_nx+1][1]);
	m_P[m_nx+1][m_ny+1] = 0.5 * (m_P[m_nx][m_ny+1] + m_P[m_nx+1][m_ny]);
}

void IncompressibleKernel::defineBody(std::vector<std::vector<bool>> v_body)
{
	m_body = v_body;
}

void IncompressibleKernel::loadFieldsAt(size_t t)
{
	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			if (m_body[i][j])
			{
				m_U[i+1][j+1] = 0.;
				U_tilde[i+1][j+1] = 0.;
			}
			else
			{
				m_U[i+1][j+1] = m_data.getXVelocityAt(t, i, j);
				U_tilde[i+1][j+1] = m_U[i+1][j+1] + m_dt * (diffusion_u(i+1, j+1) - advection_u(i+1, j+1) - pressure_u(i+1, j+1));
			}
			vec_X[i + m_nx*j] = m_U[i+1][j+1];
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			if (m_body[i][j])
			{
				m_V[i+1][j+1] = 0.;
				V_tilde[i+1][j+1] = 0.;
			}
			else
			{
				m_V[i+1][j+1] = m_data.getYVelocityAt(t, i, j);
				V_tilde[i+1][j+1] = m_V[i+1][j+1] + m_dt * (diffusion_v(i+1, j+1) - advection_v(i+1, j+1) - pressure_v(i+1, j+1));
			}
			vec_X[m_nb_pts + i + m_nx*j] = m_V[i+1][j+1];
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			m_P[i+1][j+1] = m_data.getPressureAt(t, i, j);
			vec_X[2*m_nb_pts + i + m_nx*j] = m_P[i+1][j+1];
		}
	}
}

void IncompressibleKernel::saveFieldsAt(size_t t) const
{
	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			m_data.setXVelocityAt(t, i, j, vec_X[i + m_nx*j]);
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			m_data.setYVelocityAt(t, i, j, vec_X[m_nb_pts + i + m_nx*j]);
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			m_data.setPressureAt(t, i, j, vec_X[2*m_nb_pts + i + m_nx*j]);
		}
	}
}

/*
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
*/

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
	return m_U[i][j] * dudx(i,j) + m_V[i][j] * dudy(i,j);
}

double IncompressibleKernel::advection_v(size_t i, size_t j) const
{
	return m_U[i][j] * dvdx(i,j) + m_V[i][j] * dvdy(i,j);
}

double IncompressibleKernel::diffusion_u(size_t i, size_t j) const
{
	return m_nu * laplacian_u(i,j);
}

double IncompressibleKernel::diffusion_v(size_t i, size_t j) const
{
	return m_nu * laplacian_v(i,j);
}

double IncompressibleKernel::pressure_u(size_t i, size_t j) const
{
	return (m_P[i+1][j] - m_P[i-1][j])/(2*m_dx*m_rho);
}

double IncompressibleKernel::pressure_v(size_t i, size_t j) const
{
	return (m_P[i][j+1] - m_P[i][j-1])/(2*m_dy*m_rho);
}
/*
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
*/
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
/*
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
*/

void IncompressibleKernel::computeMatrixA()
{
	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			diag_beta1[i + m_nx*j] = -m_theta*(m_nu_x - 0.5*m_dtdx*U_tilde[i+1][j+1]);
			diag_beta2[i + m_nx*j] = -m_theta*(m_nu_x + 0.5*m_dtdx*U_tilde[i+1][j+1]);
			diag_gamma1[i + m_nx*j] = -m_theta*(m_nu_y - 0.5*m_dtdy*V_tilde[i+1][j+1]);
			diag_gamma2[i + m_nx*j] = -m_theta*(m_nu_y + 0.5*m_dtdy*V_tilde[i+1][j+1]);

			diag_Au[i + m_nx*j] = -2.*m_rho*m_theta*U_tilde[i+1][j+1];
			diag_B1u[i + m_nx*j] = m_rho*m_theta*U_tilde[i+2][j+1];
			diag_B2u[i + m_nx*j] = m_rho*m_theta*U_tilde[i][j+1];
			diag_D1[i + m_nx*j] = 0.25*m_dxdy*m_rho*m_theta*V_tilde[i+2][j+2];
			diag_D2[i + m_nx*j] = -0.25*m_dxdy*m_rho*m_theta*V_tilde[i+2][j];
			diag_D3[i + m_nx*j] = -0.25*m_dxdy*m_rho*m_theta*V_tilde[i][j+2];
			diag_D4[i + m_nx*j] = 0.25*m_dxdy*m_rho*m_theta*V_tilde[i][j];

			diag_Av[i + m_nx*j] = -2.*m_rho*m_theta*m_dxdy*m_dxdy*V_tilde[i+1][j+1];
			diag_C1v[i + m_nx*j] = m_rho*m_theta*m_dxdy*m_dxdy*V_tilde[i+1][j+2];
			diag_C2v[i + m_nx*j] = m_rho*m_theta*m_dxdy*m_dxdy*V_tilde[i+1][j];
			diag_E1[i + m_nx*j] = 0.25*m_dxdy*m_rho*m_theta*U_tilde[i+2][j+2];
			diag_E2[i + m_nx*j] = -0.25*m_dxdy*m_rho*m_theta*U_tilde[i+2][j];
			diag_E3[i + m_nx*j] = -0.25*m_dxdy*m_rho*m_theta*U_tilde[i][j+2];
			diag_E4[i + m_nx*j] = 0.25*m_dxdy*m_rho*m_theta*U_tilde[i][j];
		}
	}
}

void IncompressibleKernel::computeVectorB()
{
	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			vec_B[i + m_nx*j] = m_U[i+1][j+1] + (1.-m_theta)*m_dt * (diffusion_u(i+1,j+1) - advection_u(i+1,j+1));
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			vec_B[m_nb_pts + i + m_nx*j] = m_V[i+1][j+1] + (1.-m_theta)*m_dt * (diffusion_v(i+1,j+1) - advection_v(i+1,j+1));
		}
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			vec_B[2*m_nb_pts + i + m_nx*j] = m_rho*(m_theta-1.) * (d2u2dx2(i+1,j+1) + 2.*d2uvdxdy(i+1,j+1) + d2v2dy2(i+1,j+1));
		}
	}
}

double IncompressibleKernel::safe_access_VectorX(size_t block, size_t i, size_t j) const
{
	if (i == -1)
	{
		if (j == -1)
		{
			switch (block)
			{
				case 0: return m_U[0][0];
				case 1: return m_V[0][0];
				case 2: return m_P[0][0];
			}
		}
		else if (j == m_ny)
		{
			switch (block)
			{
				case 0: return m_U[0][m_ny+1];
				case 1: return m_V[0][m_ny+1];
				case 2: return m_P[0][m_ny+1];
			}
		}
		else
		{
			switch (block)
			{
				case 0: return m_U[0][j+1];
				case 1: return m_V[0][j+1];
				case 2: return m_P[0][j+1];
			}
		}
	}

	if (i == m_nx)
	{
		if (j == -1)
		{
			switch (block)
			{
				case 0: return m_U[m_nx+1][0];
				case 1: return m_V[m_nx+1][0];
				case 2: return m_P[m_nx+1][0];
			}
		}
		else if (j == m_ny)
		{
			switch (block)
			{
				case 0: return m_U[m_nx+1][m_ny+1];
				case 1: return m_V[m_nx+1][m_ny+1];
				case 2: return m_P[m_nx+1][m_ny+1];
			}
		}
		else
		{
			switch (block)
			{
				case 0: return m_U[m_nx+1][j+1];
				case 1: return m_V[m_nx+1][j+1];
				case 2: return m_P[m_nx+1][j+1];
			}
		}
	}

	if (j == -1)
	{
		switch (block)
		{
			case 0: return m_U[i+1][0];
			case 1: return m_V[i+1][0];
			case 2: return m_P[i+1][0];
		}
	}

	if (j == m_ny)
	{
		switch (block)
		{
			case 0: return m_U[i+1][m_ny+1];
			case 1: return m_V[i+1][m_ny+1];
			case 2: return m_P[i+1][m_ny+1];
		}
	}

	return vec_X[block*m_nb_pts + i + m_nx*j];
}

double IncompressibleKernel::GS_iteration_U()
{
	double mean_diff(0.);
	double new_Xk;

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			new_Xk = vec_B[i + m_nx*j];

			new_Xk -= diag_gamma2[i + m_nx*j] * safe_access_VectorX(0, i, j-1);
			new_Xk -= diag_beta2[i + m_nx*j] * safe_access_VectorX(0, i-1, j);
			new_Xk -= diag_beta1[i + m_nx*j] * safe_access_VectorX(0, i+1, j);
			new_Xk -= diag_gamma1[i + m_nx*j] * safe_access_VectorX(0, i, j+1);
			new_Xk -= diag_delta2u[i + m_nx*j] * safe_access_VectorX(2, i-1, j);
			new_Xk -= diag_delta1u[i + m_nx*j] * safe_access_VectorX(2, i+1, j);

			new_Xk /= diag_alpha[i + m_nx*j];

			mean_diff += new_Xk - vec_X[i + m_nx*j];
			vec_X[i + m_nx*j] = new_Xk;
		}
	}

	return mean_diff / m_nb_pts;
}

double IncompressibleKernel::GS_iteration_V()
{
	double mean_diff(0.);
	double new_Xk;

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			new_Xk = vec_B[m_nb_pts + i + m_nx*j];

			new_Xk -= diag_gamma2[i + m_nx*j] * safe_access_VectorX(1, i, j-1);
			new_Xk -= diag_beta2[i + m_nx*j] * safe_access_VectorX(1, i-1, j);
			new_Xk -= diag_beta1[i + m_nx*j] * safe_access_VectorX(1, i+1, j);
			new_Xk -= diag_gamma1[i + m_nx*j] * safe_access_VectorX(1, i, j+1);
			new_Xk -= diag_delta2v[i + m_nx*j] * safe_access_VectorX(2, i, j-1);
			new_Xk -= diag_delta1v[i + m_nx*j] * safe_access_VectorX(2, i, j+1);
			
			new_Xk /= diag_alpha[i + m_nx*j];

			mean_diff += new_Xk - vec_X[m_nb_pts + i + m_nx*j];
			vec_X[m_nb_pts + i + m_nx*j] = new_Xk;
		}
	}

	return mean_diff / m_nb_pts;
}

double IncompressibleKernel::GS_iteration_P()
{
	double mean_diff(0.);
	double new_Xk;

	for (size_t j = 0; j < m_ny; j++)
	{
		for (size_t i = 0; i < m_nx; i++)
		{
			new_Xk = vec_B[2*m_nb_pts + i + m_nx*j];

			new_Xk -= diag_D4[i + m_nx*j] * safe_access_VectorX(0, i-1, j-1);
			new_Xk -= diag_D2[i + m_nx*j] * safe_access_VectorX(0, i+1, j-1);
			new_Xk -= diag_B2u[i + m_nx*j] * safe_access_VectorX(0, i-1, j);
			new_Xk -= diag_Au[i + m_nx*j] * safe_access_VectorX(0, i, j);
			new_Xk -= diag_B1u[i + m_nx*j] * safe_access_VectorX(0, i+1, j);
			new_Xk -= diag_D3[i + m_nx*j] * safe_access_VectorX(0, i-1, j+1);
			new_Xk -= diag_D1[i + m_nx*j] * safe_access_VectorX(0, i+1, j+1);

			new_Xk -= diag_E4[i + m_nx*j] * safe_access_VectorX(1, i-1, j-1);
			new_Xk -= diag_C2v[i + m_nx*j] * safe_access_VectorX(1, i, j-1);
			new_Xk -= diag_E2[i + m_nx*j] * safe_access_VectorX(1, i+1, j-1);
			new_Xk -= diag_Av[i + m_nx*j] * safe_access_VectorX(1, i, j);
			new_Xk -= diag_E3[i + m_nx*j] * safe_access_VectorX(1, i-1, j+1);
			new_Xk -= diag_C1v[i + m_nx*j] * safe_access_VectorX(1, i, j+1);
			new_Xk -= diag_E1[i + m_nx*j] * safe_access_VectorX(1, i+1, j+1);

			new_Xk -= diag_C2[i + m_nx*j] * safe_access_VectorX(2, i, j-1);
			new_Xk -= diag_B2[i + m_nx*j] * safe_access_VectorX(2, i-1, j);
			new_Xk -= diag_B1[i + m_nx*j] * safe_access_VectorX(2, i+1, j);
			new_Xk -= diag_C1[i + m_nx*j] * safe_access_VectorX(2, i, j+1);

			new_Xk /= diag_A[i + m_nx*j];

			mean_diff += new_Xk - vec_X[2*m_nb_pts + i + m_nx*j];
			vec_X[2*m_nb_pts + i + m_nx*j] = new_Xk;
		}
	}

	return mean_diff / m_nb_pts;
}

void IncompressibleKernel::solve_GaussSeidel()
{
	double mean_diff_u, mean_diff_v, mean_diff_p;
	bool converged;

	do
	{
		mean_diff_u = GS_iteration_U();
		mean_diff_v = GS_iteration_V();
		mean_diff_p = GS_iteration_P();

		converged = (mean_diff_u < m_accuracy_u) && (mean_diff_v < m_accuracy_v) && (mean_diff_p < m_accuracy_p);
	}
	while (!converged);
}

void IncompressibleKernel::simulate_theta()
{
	loadFieldsAt(0);
	PoissonSolver pressure_solver(m_Lx, m_Ly, m_nx, m_ny);
	std::cout << "Computing the initial pressure field..." << std::endl;
	computePressureField(pressure_solver, 0);

	std::cout << "Simulating the incompressible flow..." << std::endl;

	for (size_t t = 0; t < m_nb_steps; t++)
	{
		loadFieldsAt(t);
		printSimulationProgression(t);

		computeMatrixA();
		computeVectorB();
		solve_GaussSeidel();

		saveFieldsAt(t+1);
	}

	std::cout << "  --->  OK" << std::endl;

	if (m_clean_pressure)
	{
		std::cout << "Recomputing the pressure..." << std::endl;
		
		for (size_t t = 0; t < m_nb_steps; t++)
		{
			printSimulationProgression(t);
			computePressureField(pressure_solver, t);
		}
		
		std::cout << "  --->  OK" << std::endl;
	}
}