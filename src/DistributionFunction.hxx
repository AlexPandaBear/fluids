#pragma once

#include <memory>

/*! @class DistributionFunction
 *
 *  @brief A LBM distribution function
 *	
 *	This class represents the state of a discrete distribution function.
 *  The purpose of an instance of this class is to overload the () operator, in order to ease the writing of the LBM.
 *	This is a 3 dimensional object : two dimensions for space and one for velocities.
 *	It is therefore only designed to represent 2D instantaneous distributions on cartesian grids, and can only be interpreted with a velocity scheme.
 *  The values are stored in a dynamically allocated array of double
 */
class DistributionFunction
{
private:
	size_t m_nx, m_ny, m_nv;
	std::unique_ptr<double[]> ptr_f;

public:
	/*! @brief The constructor of the class
	 *
	 *  This constructor dynamically allocates the memory needed to store a field, which dimensions are defined by the parameters.
	 *
	 *	@param nx The number of points along the x-axis
	 *
	 *	@param ny The number of points along the y-axis
	 *
	 *	@param nv The number of velocities if the velocity scheme
	 *
	 *  @warning This constructor only allocates the memory, but does not initialize the values of the field.
	 */
	DistributionFunction(size_t nx, size_t ny, size_t nv);

	//DistributionFunction(std::string file_name);
	
	/*! @brief The destructor of the class
	 *
	 */
	~DistributionFunction();

	double& operator()(size_t i, size_t j, size_t k);

	//void save(std::string file_name) const;
};


double equilibrium(..);