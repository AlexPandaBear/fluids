#include <iostream>
#include "DistributionFunction.hxx"

int main(int argc, char const *argv[])
{
	size_t nx(2), ny(3), nv(5);

	//Test constructor
	DistributionFunction f(nx, ny, nv);

	//Test velocity scheme
	f.set_velocity_scheme({0, 1, 0, -1, 0}, {0, 0, 1, 0, -1}, {0.5, 0.125, 0.125, 0.125, 0.125});
	f.display_velocity_scheme();

	//Test operator()
	for (size_t i = 0; i < nx; i++)
	{
		for (size_t j = 0; j < ny; j++)
		{
			for (size_t k = 0; k < nv; k++)
			{
				f(i, j, k) = i + j + k;
			}
		}
	}

	for (size_t k = 0; k < nv; k++)
	{
		for (size_t i = 0; i < nx; i++)
		{
			for (size_t j = 0; j < ny; j++)
			{
				std::cout << f(i, j, k) << " ";
			}

			std::cout << std::endl;
		}

		std::cout << std::endl;
	}

	return 0;
}