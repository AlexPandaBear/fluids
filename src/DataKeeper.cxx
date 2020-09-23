#include "DataKeeper.hxx"

DataKeeper::DataKeeper() :
	m_nb_steps(0),
	m_nx(0),
	m_ny(0),
	m_coef_1(0.),
	m_coef_2(0.),
	ptr_T(NULL),
	ptr_P(NULL),
	ptr_U(NULL),
	ptr_V(NULL) {}

DataKeeper::~DataKeeper() {}

void DataKeeper::reset_size(size_t nb_steps, size_t nx, size_t ny)
{
	m_nb_steps = nb_steps;
	m_nx = nx;
	m_ny = ny;

	m_coef_1 = nx*ny;
	m_coef_2 = ny;

	ptr_T.reset(new double[(nb_steps+1)*nx*ny]);
	ptr_P.reset(new double[(nb_steps+1)*nx*ny]);
	ptr_U.reset(new double[(nb_steps+1)*nx*ny]);
	ptr_V.reset(new double[(nb_steps+1)*nx*ny]);
}

double DataKeeper::getTemperatureAt(size_t t, size_t i, size_t j) const
{
	return ptr_T[m_coef_1*t + m_coef_2*i + j];
}

double DataKeeper::getPressureAt(size_t t, size_t i, size_t j) const
{
	return ptr_P[m_coef_1*t + m_coef_2*i + j];
}

double DataKeeper::getXVelocityAt(size_t t, size_t i, size_t j) const
{
	return ptr_U[m_coef_1*t + m_coef_2*i + j];
}

double DataKeeper::getYVelocityAt(size_t t, size_t i, size_t j) const
{
	return ptr_V[m_coef_1*t + m_coef_2*i + j];
}

void DataKeeper::setTemperatureAt(size_t t, size_t i, size_t j, double T)
{
	ptr_T[m_coef_1*t + m_coef_2*i + j] = T;
}

void DataKeeper::setPressureAt(size_t t, size_t i, size_t j, double P)
{
	ptr_P[m_coef_1*t + m_coef_2*i + j] = P;
}

void DataKeeper::setXVelocityAt(size_t t, size_t i, size_t j, double U)
{
	ptr_U[m_coef_1*t + m_coef_2*i + j] = U;
}

void DataKeeper::setYVelocityAt(size_t t, size_t i, size_t j, double V)
{
	ptr_V[m_coef_1*t + m_coef_2*i + j] = V;
}