#include "DataProcessor.hxx"

DataProcessor::DataProcessor(DataKeeper& data) :
	m_ready(false),
	m_data(data) {}

DataProcessor::~DataProcessor() {}

bool DataProcessor::is_ready() const
{
	return m_ready;
}

void DataProcessor::define_parameters(DataKeeper &data, size_t nx, size_t ny, size_t nb_steps, double dx, double dy, double dt)
{
	m_nx = nx;
	m_ny = ny;
	m_nb_steps = nb_steps;
	m_dx = dx;
	m_dy = dy;
	m_dt = dt;
	m_ready = true;
}

double DataProcessor::dudx(size_t t, size_t i, size_t j) const
{
	return (m_data.getXVelocityAt(t, i+1, j) - m_data.getXVelocityAt(t, i-1, j)) / 2.*m_dx;
}

double DataProcessor::dvdy(size_t t, size_t i, size_t j) const
{
	return (m_data.getYVelocityAt(t, i, j+1) - m_data.getYVelocityAt(t, i, j-1)) / 2.*m_dy;
}

std::vector<std::vector<double>> DataProcessor::norm_V_field(size_t t) const
{
	std::vector<std::vector<double>> norm_V(m_nx, std::vector<double>(m_ny, 0.));

	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			norm_V[i][j] = sqrt( pow(m_data.getXVelocityAt(t, i, j), 2) + pow(m_data.getYVelocityAt(t, i, j), 2) );
		}
	}

	return norm_V;
}

std::vector<std::vector<double>> DataProcessor::div_V_field(size_t t) const
{
	std::vector<std::vector<double>> div_V(m_nx, std::vector<double>(m_ny, 0.));

	for (size_t i = 1; i < m_nx-1; i++)
	{
		for (size_t j = 1; j < m_ny-1; j++)
		{
			div_V[i][j] = dudx(t, i, j) + dvdy(t, i, j);
		}
	}

	return div_V;
}