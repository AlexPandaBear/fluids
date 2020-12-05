#include "Field.hxx"

Field::Field(std::vector<size_t> dimensions) :
	m_nb_dim(dimensions.size()),
	m_dimensions(dimensions),
	ptr_data(NULL)
{
	size_t size = 1;

	for (size_t d = 0; d < m_nb_dim; d++)
	{
		size *= m_dimensions[d];
	}

	ptr_data.reset(new double[size]);
}

Field::~Field() {}

double& Field::operator()(std::vector<size_t> index)
{
	size_t raw_index = index[0];

	for (size_t d = 1; d < m_nb_dim; d++)
	{
		raw_index *= m_dimensions[d];
		raw_index += index[d];
	}

	return ptr_data[raw_index];
}