#include "LBMKernel.hxx"

int main(int argc, char const *argv[])
{
	LBMKernel lbm;

	lbm.build_space(0., 1., 3, 0., 1., 3);
	lbm.build_time(0., 1., 3);

	std::vector<std::pair<size_t, size_t>> ij = { {0,0}, {0,1}, {0,2}, {1,0}, {1,2}, {2,0}, {2,1}, {2,2} };
	std::vector<std::vector<size_t>> dir = { {2,3,4,5,6}, {2,3,4}, {8,1,2,3,4}, {4,5,6}, {8,1,2}, {4,5,6,7,8}, {6,7,8}, {6,7,8,1,2} };

	D2Q9BounceBack walls(ij, dir);

	BGK c(1., 1.4);
	Stream2D s({walls});
	Field f;
	lbm.simulate(c, s, f);

	return 0;
}