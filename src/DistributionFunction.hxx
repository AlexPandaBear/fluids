#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>

/*! @class DistributionFunction
 *
 *  @brief A LBM distribution function
 *	
 *	This class represents the state of a discrete distribution function.
 *  The purpose of an instance of this class is to overload the () operator, in order to ease the writing of the LBM.
 *	This is a 3 dimensional object : two dimensions for space and one for velocities.
 *	It is therefore only designed to represent 2D instantaneous distributions on cartesian grids, and can only be interpreted with a velocity scheme.
 *  The values are stored in a dynamically allocated array of double.
 */
class DistributionFunction
{
private:
	size_t m_nb_dim;
	std::vector<size_t> m_sizes;
	
	size_t m_nv;
	std::vector<std::vector<int>> m_e;
	std::vector<double> m_e2, m_w;

	std::unique_ptr<double[]> ptr_f;

public:
	/*! @brief The constructor of the class
	 *
	 *  This constructor dynamically allocates the memory needed to store a field, which dimensions are defined by the parameters.
	 *
	 *	@param sizes The number of points along each axis
	 *
	 *	@param vs The velocity scheme to use
	 *
	 *  @warning This constructor only allocates the memory, but does not initialize the values of the field.
	 *
	 *	@warning This constructor only allocates the right amount of memory given the size of the velocity scheme,
	 *	but the velocity scheme should be defined with the set_velocity_scheme method before use.
	 */
	DistributionFunction(std::vector<size_t> sizes, std::string vs);

	//DistributionFunction(std::string file_name);
	
	/*! @brief The destructor of the class
	 *
	 */
	~DistributionFunction();

	/*! @brief A method to define the velocity scheme
	 *
	 *	@param elem_vel A std::vector whose elements descibe the velocities of the scheme
	 *
	 *	@param weights A std::vector defining the weight of each one of the velocities for gaussian quadrature
	 */
	void set_velocity_scheme(std::vector<std::vector<int>> elem_vel, std::vector<double> weights);

	//void display_velocity_scheme() const;

	/*! @brief A method to get the dimensions of the mesh
	 *
	 *	@returns A std::vector of the number of points in each direction
	 */
	std::vector<size_t> get_sizes() const;

	/*! @brief A method to get one velocity of the scheme
	 *
	 *	@param k The id of the velocity in the scheme
	 *
	 *	@returns A std::vector of the components of the velocity
	 */
	std::vector<int> get_e(size_t k) const;

	/*! @brief A method to get the gaussian quadrature weight of one velocity of the scheme
	 *
	 *	@param k The id of the velocity in the scheme
	 *
	 *	@returns The weight associated to the velocity
	 */
	double get_w(size_t k) const;

	/*! @brief A method to get the state of the distribution function at a given point in space
	 *
	 *	@param pt The index of the point on each axis
	 *
	 *	@returns The value associated to each elementary velocity at this point
	 */
	std::vector<double> get_f(std::vector<size_t> pt) const;

	/*! @brief A method computing the density
	 *
	 *	This method computes the density of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of the distribution function.
	 *
	 *	@param pt The index of the point on each axis
	 *
	 *	@returns The fluid mass by volume unit
	 */
	double compute_rho(std::vector<size_t> pt) const;

	/*! @brief A method computing the momentum along the x-axis
	 *
	 *	This method computes the projection on a specific axis of the momentum of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of the distribution function multiplied by the microscopic velocities.
	 *
	 *	@param pt The index of the point on each axis
	 *
	 *	@param axis The id of the axis to project on
	 *
	 *	@returns The fluid momentum projected on the chosen axis by volume unit
	 */
	double compute_rho_ui(std::vector<size_t> pt, size_t axis) const;

	/*! @brief A method computing the energy
	 *
	 *	This method computes the energy of the fluid at a specific point in space.
	 *	It implements the gaussian quadrature over the discrete velocity space of half the distribution function multiplied by the square of the microscopic velocities.
	 *
	 *	@param pt The index of the point on each axis
	 *
	 *	@returns The fluid energy by volume unit
	 */
	double compute_rho_e(std::vector<size_t> pt) const;

	/*! @brief A method computing all the macroscopic variables
	 *
	 *	This method computes the density, all the components of the velocity and the energy per mass unit at a specific point in space.
	 *	It implements the gaussian quadratures of the distribution function's moments.
	 *
	 *	@param pt The index of the point on each axis
	 *
	 *	@returns A std::vector all the macroscopic variables listed above, in that specific order.
	 */
	std::vector<double> compute_macro_variables(std::vector<size_t> pt) const;

	/*! @brief A method to set the local values of the function
	 *
	 *	This method sets the different values of the distribution function at a given point in space.
	 *	It uses the macroscopic flow variables and assumes local equilibrium to compute the values.
	 *
	 *	@param pt The index of the point on each axis
	 *
	 *	@param rho The density at this point
	 *
	 *	@param u The velocity at this point
	 *
	 *	@param e The energy at this point
	 *
	 *	@param gamma The adiabatic coefficient of the fluid
	 *
	 *	@warning Calling this method will overwrite the previous values of the distribution function at the point considered.
	 */
	void set_equilibrium(std::vector<size_t> pt, double rho, std::vector<double> u, double e, double gamma);

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