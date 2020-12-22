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

void Field::reshape(std::vector<size_t> dimensions)
{
	m_nb_dim = dimensions.size();
	m_dimensions = dimensions;

	size_t size = 1;

	for (size_t d = 0; d < m_nb_dim; d++)
	{
		size *= m_dimensions[d];
	}

	ptr_data.reset(new double[size]);
}

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

double* Field::get_ptr() const
{
	return ptr_data.get();
}

size_t Field::get_element_size() const
{
	return sizeof(double);
}

size_t Field::get_nb_dimensions() const
{
	return m_nb_dim;
}

std::vector<size_t> Field::get_dimensions_sizes() const
{
	return m_dimensions;
}

std::vector<size_t> Field::get_adjusted_element_sizes() const
{
	std::vector<size_t> res(m_nb_dim, get_element_size());

	for (size_t d = 0; d < m_nb_dim; d++)
	{
		for (size_t i = 0; i < d; i++)
		{
			res[i] *= m_dimensions[d];
		}
	}

	return res;
}