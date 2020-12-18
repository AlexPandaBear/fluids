#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <string>

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

	std::vector<int> m_ex, m_ey;
	std::vector<double> m_e2, m_w;

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
	 *
	 *	@warning This constructor only allocates the right amount of memory given the size of the velocity scheme,
	 *	but the velocity scheme should be defined with the set_velocity_scheme method before use.
	 */
	DistributionFunction(size_t nx = 0, size_t ny = 0, size_t nv = 0);

	DistributionFunction(size_t nx, size_t ny, size_t nv, std::string vs);

	//DistributionFunction(std::string file_name);
	
	/*! @brief The destructor of the class
	 *
	 */
	~DistributionFunction();

	/*! @brief A method to define the velocity scheme
	 *
	 *	@param ex A std::vector describing the x component of each one of the velocities
	 *
	 *	@param ey A std::vector describing the y component of each one of the velocities
	 *
	 *	@param weights A std::vector defining the weight of each one of the velocities for gaussian quadrature
	 */
	void set_velocity_scheme(std::vector<int> ex, std::vector<int> ey, std::vector<double> weights);

	/*! @brief A method displaying the velocity scheme
	 */
	void display_velocity_scheme() const;

	/*! @brief A method to get the number of points in the x direction
	 */
	size_t get_size_x() const;
	
	/*! @brief A method to get the number of points in the y direction
	 */
	size_t get_size_y() const;

	/*! @brief A method to get the number of discrete velocities
	 */
	size_t get_size_v() const;

	/*! @brief A method to get the x component of one velocity of the scheme
	 *
	 *	@param k The id of the velocity in the scheme
	 */
	int get_ex(size_t k) const;
	
	/*! @brief A method to get the y component of one velocity of the scheme
	 *
	 *	@param k The id of the velocity in the scheme
	 */
	int get_ey(size_t k) const;

	/*! @brief A method to get the gaussian quadrature weight of one velocity of the scheme
	 *
	 *	@param k The id of the velocity in the scheme
	 */
	double get_w(size_t k) const;

	/*! @brief A method computing the density
	 *
	 *	This method computes the density of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of the distribution function.
	 *
	 *	@param i The index along the x-axis of the point
	 *
	 *	@param j The index along the y-axis of the point
	 *
	 *	@returns The fluid mass by volume unit
	 */
	double compute_rho(size_t i, size_t j) const;

	/*! @brief A method computing the momentum along the x-axis
	 *
	 *	This method computes the projection on the x-axis of the momentum of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of the distribution function multiplied by the microscopic velocities.
	 *
	 *	@param i The index along the x-axis of the point
	 *
	 *	@param j The index along the y-axis of the point
	 *
	 *	@returns The fluid x-momentum by volume unit
	 */
	double compute_rho_ux(size_t i, size_t j) const;

	/*! @brief A method computing the momentum along the y-axis
	 *
	 *	This method computes the projection on the y-axis of the momentum of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of the distribution function multiplied by the microscopic velocities.
	 *
	 *	@param i The index along the x-axis of the point
	 *
	 *	@param j The index along the y-axis of the point
	 *
	 *	@returns The fluid y-momentum by volume unit
	 */
	double compute_rho_uy(size_t i, size_t j) const;

	/*! @brief A method computing the energy
	 *
	 *	This method computes the energy of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of half the distribution function multiplied by the square of the microscopic velocities.
	 *
	 *	@param i The index along the x-axis of the point
	 *
	 *	@param j The index along the y-axis of the point
	 *
	 *	@returns The fluid energy by volume unit
	 */
	double compute_rho_e(size_t i, size_t j) const;

	/*! @brief An operator serving as a getter/setter for the values of the distribution function
	 *
	 *	This operator is similar to the [] operator for arrays or vectors.
	 *	It allows to get or to set one specific value of the function.
	 *	
	 *	@param i The index along the x-axis of the point
	 *
	 *	@param j The index along the y-axis of the point
	 *
	 *	@param k The velocity number at the point
	 *
	 *	@returns The reference to the value
	 */
	double& operator()(size_t i, size_t j, size_t k);

	//void save(std::string file_name) const;
};