#include <iostream>
#include "Field.hxx"

int main(int argc, char const *argv[])
{
	//Test constructor
	Field f({5, 3, 2});
	//std::cout << "coucou" << std::endl;

	//Test operator()
	for (size_t i = 0; i < 5; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			for (size_t k = 0; k < 2; k++)
			{
				f({i, j, k}) = i+j;
			}
		}
	}


	for (size_t i = 0; i < 5; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			for (size_t k = 0; k < 2; k++)
			{
				std::cout << f({i, j, k}) << " ";
			}
		}
	}

	std::cout << std::endl;

	return 0;
}