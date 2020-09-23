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

void DataKeeper::saveData(std::string file_name) const
{
	std::ofstream outfile;
	outfile.open(file_name);

	if (!outfile.is_open())
	{
		throw std::invalid_argument("Unable to save data in " + file_name);
	}

	outfile << m_nb_steps << "\n";
	outfile << m_nx << "\n";
	outfile << m_ny << "\n";

	size_t index_max((m_nb_steps+1)*m_nx*m_ny);

	for (size_t index = 0; index < index_max; index++)
	{
		outfile << ptr_T[index] << "\n";
	}

	for (size_t index = 0; index < index_max; index++)
	{
		outfile << ptr_P[index] << "\n";
	}

	for (size_t index = 0; index < index_max; index++)
	{
		outfile << ptr_U[index] << "\n";
	}

	for (size_t index = 0; index < index_max; index++)
	{
		outfile << ptr_V[index] << "\n";
	}

	outfile.close();
}

void DataKeeper::loadData(std::string file_name)
{
	std::ifstream file;
	file.open(file_name);

	if (!file.is_open())
	{
		throw std::invalid_argument("Unable to read data in " + file_name);
	}

	std::string data_str;

	getline(file, data_str);
	m_nb_steps = atoi(data_str.c_str());

	getline(file, data_str);
	m_nx = atoi(data_str.c_str());

	getline(file, data_str);
	m_ny = atoi(data_str.c_str());

	reset_size(m_nb_steps, m_nx, m_ny);
	size_t index_max((m_nb_steps+1)*m_nx*m_ny);

	for (size_t index = 0; index < index_max; index++)
	{
		getline(file, data_str);
		ptr_T[index] = atof(data_str.c_str());
	}

	for (size_t index = 0; index < index_max; index++)
	{
		getline(file, data_str);
		ptr_P[index] = atof(data_str.c_str());
	}

	for (size_t index = 0; index < index_max; index++)
	{
		getline(file, data_str);
		ptr_U[index] = atof(data_str.c_str());
	}

	for (size_t index = 0; index < index_max; index++)
	{
		getline(file, data_str);
		ptr_V[index] = atof(data_str.c_str());
	}

	file.close();
}