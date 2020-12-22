#include "LBMKernel.hxx"

LBMKernel::LBMKernel() {}

LBMKernel::~LBMKernel() {}

void LBMKernel::build_space(double x_min, double x_max, size_t nx, double y_min, double y_max, size_t ny)
{
	m_x = std::vector<double>(nx, 0.);
	m_y = std::vector<double>(ny, 0.);

	double dx((x_max - x_min) / (nx-1));
	double dy((y_max - y_min) / (ny-1));

	for (size_t i = 0; i < nx; i++)
	{
		m_x[i] = x_min + i*dx;
	}

	for (size_t j = 0; j < nx; j++)
	{
		m_y[j] = y_min + j*dy;
	}
}

void LBMKernel::build_time(double t_start, double t_end, size_t n)
{
	m_time = std::vector<double>(n+1, 0.);

	double dt((t_end - t_start) / n);

	for (size_t t = 0; t < n+1; t++)
	{
		m_time[t] = t_start + t*dt;
	}
}

void LBMKernel::set_initial_state(Field& rho, Field& U, Field& E, double gamma)
{
	size_t nx(m_f.get_size_x());
	size_t ny(m_f.get_size_y());

	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			m_f.set_equilibrium(i, j, rho({i,j}), U({i,j,0}), U({i,j,1}), E({i,j}), gamma);
		}
	}
}

void LBMKernel::simulate(CollisionOperator& collision, StreamingOperator& streaming, Field& f_out)
{
	size_t nb_steps(m_time.size() - 1);
	size_t nx(m_x.size());
	size_t ny(m_y.size());

	f_out.reshape({nb_steps, nx, ny, 9});
	streaming.initialize(m_f);

	for (size_t t = 0; t < nb_steps; t++)
	{
		collision(m_f, m_time[t+1] - m_time[t]);
		streaming(m_f);
	}
}