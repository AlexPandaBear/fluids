#include "PoissonSolver.hxx"

PoissonSolver::PoissonSolver(double Lx, double Ly, size_t nx, size_t ny) :
	m_Lx(Lx),
	m_Ly(Ly),
	m_nx(nx),
	m_ny(ny),
	m_nb_eq(nx*ny),
	m_dx(Lx/(nx-1)),
	m_dy(Ly/(ny-1)),
	m_coef_1(m_dx*m_dx/(m_dy*m_dy)),
	m_coef_2(-2.*(1+m_coef_1))
{
	std::cout << "Generating LU decomposition for the Poisson problem...     ";
	generate_LU();
	std::cout << "Done." << std::endl;
}

PoissonSolver::~PoissonSolver() {}

void PoissonSolver::definePoissonFunction(std::vector<double> v_values)
{
	if (v_values.size() != m_nx*m_ny)
	{
		throw std::invalid_argument("Size mismatch");
	}

	v_f = std::vector<double>(m_nx*m_ny, 0.);

	for (size_t i = 0; i < m_nb_eq; i++)
	{
		v_f[i] = v_values[i];
	}
}

void PoissonSolver::defineBoundaryConditions(double uniform_value)
{
	for (size_t i = 0; i < m_nx; i++)
	{
		v_f[i] -= m_coef_1*uniform_value;
		v_f[m_nb_eq-1 - i] -= m_coef_1*uniform_value;
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		v_f[m_nx*j] -= uniform_value;
		v_f[m_nx*(j+1) - 1] -= uniform_value;
	}
}

void PoissonSolver::defineBoundaryConditions(std::vector<double> v_up, std::vector<double> v_down, std::vector<double> v_left, std::vector<double> v_right)
{
	if (v_up.size() != m_nx || v_down.size() != m_nx || v_left.size() != m_ny || v_right.size() != m_ny)
	{
		throw std::invalid_argument("Size mismatch");
	}

	for (size_t i = 0; i < m_nx; i++)
	{
		v_f[i] -= m_coef_1*v_up[i];
		v_f[m_nb_eq-1 - i] -= m_coef_1*v_down[i];
	}

	for (size_t j = 0; j < m_ny; j++)
	{
		v_f[m_nx*j] -= v_left[j];
		v_f[m_nx*(j+1) - 1] -= v_right[j];
	}
}

double PoissonSolver::matrix_coef(size_t i, size_t k)
{
	std::vector<size_t> v_zeros(m_ny, 0.);
	
	for (size_t i = 1; i < m_ny+1; i++)
	{
		v_zeros[i] = i*m_nx;
	}

	if (abs(i-k) == m_nx)
	{
		return 1.;
	}

	bool i_in_zeros, k_in_zeros;

	switch (abs(i-k))
	{
		case 0:
			return m_coef_2;
		
		case 1:
			i_in_zeros = (std::find(v_zeros.begin(), v_zeros.end(), i) != v_zeros.end());
			k_in_zeros = (std::find(v_zeros.begin(), v_zeros.end(), k) != v_zeros.end());
			
			if (i_in_zeros && (i > k))
			{
				return 0.;
			}
			else if (k_in_zeros && (k > i))
			{
				return 0.;
			}
			
			return m_coef_1;
		
		default:
			return 0.;
	}
}

void PoissonSolver::generate_LU()
{
	mat_upper = std::vector<std::vector<double>>(m_nb_eq, std::vector<double>(m_nb_eq, 0.));
	mat_lower = std::vector<std::vector<double>>(m_nb_eq, std::vector<double>(m_nb_eq, 0.));
	
	for (size_t i = 0; i < m_nb_eq; i++)
	{
		mat_lower[i][i] = 1.;
	}

	double lu_sum;

	for (size_t i = 0; i < m_nb_eq; i++)
	{
		for (size_t k = i; k < m_nb_eq; k++)
		{
			lu_sum = 0.;

			for (size_t j = 0; j < i; j++)
			{
				lu_sum += mat_lower[i][j] * mat_upper[j][k];
			}

			mat_upper[i][k] = matrix_coef(i,k) - lu_sum;
		}

		for (size_t k = i; k < m_nb_eq; k++)
		{
			lu_sum = 0.;

			for (size_t j = 0; j < i; j++)
			{
				lu_sum += mat_lower[k][j] * mat_upper[j][i];
			}

			mat_lower[k][i] = (matrix_coef(k,i) - lu_sum)/mat_upper[i][i];
		}
	}
}

void PoissonSolver::solve_for_y()
{
	v_y = std::vector<double>(m_nb_eq, 0.);

	double ly_sum;

	for (size_t i = 0; i < m_nb_eq; i++)
	{
		ly_sum = 0.;

		for (size_t k = 0; k < i; k++)
		{
			ly_sum += mat_lower[i][k]*v_y[k];
		}

		v_y[i] = v_f[i] - ly_sum;
	}
}

void PoissonSolver::solve_for_x()
{
	v_x = std::vector<double>(m_nb_eq, 0.);

	double ux_sum;
	size_t h;

	for (size_t i = 0; i < m_nb_eq ; i++)
	{
		h = m_nb_eq-1 - i;
		ux_sum = 0.;

		for (size_t k = h+1; k < m_nb_eq; k++)
		{
			ux_sum += mat_upper[h][k]*v_x[k];
		}

		v_x[h] = (v_y[h] - ux_sum)/mat_upper[h][h];
	}
}

std::vector<double> PoissonSolver::solve()
{
	solve_for_y();
	solve_for_x();
	return v_x;
}