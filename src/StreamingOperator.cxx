#include "StreamingOperator.hxx"

StreamingOperator::~StreamingOperator() {}

void StreamingOperator::initialize(DistributionFunction const& f) {}

void StreamingOperator::operator()(DistributionFunction& f) {}

Stream2D::Stream2D(std::vector<BoundaryCondition> BC) :
	m_nb_BC(BC.size()),
	m_BC(BC) {}

Stream2D::~Stream2D() {}

void Stream2D::initialize(DistributionFunction const& f)
{
	m_nx = f.get_size_x();
	m_ny = f.get_size_y();
	m_nv = f.get_size_v();
	m_tmp.reshape({m_nx+2, m_ny+2, m_nv});
}

void Stream2D::operator()(DistributionFunction& f)	
{
	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			for (size_t k = 0; k < m_nv; k++)
			{
				m_tmp({i+1 + f.get_ex(k), j+1 + f.get_ey(k), k}) = f(i, j, k);
			}
		}
	}

	for (size_t i = 0; i < m_nb_BC; i++)
	{
		m_BC[i].apply(m_tmp);
	}

	for (size_t i = 0; i < m_nx; i++)
	{
		for (size_t j = 0; j < m_ny; j++)
		{
			for (size_t k = 0; k < m_nv; k++)
			{
				f(i, j, k) = m_tmp({i+1, j+1, k});
			}
		}
	}
}