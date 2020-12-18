#include "BoundaryCondition.hxx"

BoundaryCondition::~BoundaryCondition() {}
	
void BoundaryCondition::apply(Field& f) const {}

D2Q9BounceBack::D2Q9BounceBack(std::vector<std::pair<size_t, size_t>> const& ij, std::vector<std::vector<size_t>> const& dir) :
	m_ij(ij),
	m_dir(dir) {}

D2Q9BounceBack::~D2Q9BounceBack() {}

void D2Q9BounceBack::apply(Field& f) const
{
	size_t i, j;
	std::vector<size_t> dir;

	for (size_t n = 0; n < m_nb_nodes; n++)
	{
		i = m_ij[n].first+1;
		j = m_ij[n].second+1;
		dir = m_dir[n];

		for (size_t d : dir)
		{
			switch (d)
			{
				case 1: f({i, j, 5}) = f({i, j+1, 1});			//	4-----3-----2
				case 2: f({i, j, 6}) = f({i-1, j+1, 2});		//	| *   |   * |
				case 3: f({i, j, 7}) = f({i-1, j, 3});			//	|   * | *   |
				case 4: f({i, j, 8}) = f({i-1, j-1, 4});		//	5-----0-----1
				case 5: f({i, j, 1}) = f({i, j-1, 5});			//	|   * | *   |
				case 6: f({i, j, 2}) = f({i+1, j-1, 6});		//	| *   |   * |
				case 7: f({i, j, 3}) = f({i+1, j, 7});			//	6-----7-----8
				case 8: f({i, j, 4}) = f({i+1, j+1, 8});
				default: throw std::logic_error("Unexpected bounce back direction");
			}
		}
	}
}