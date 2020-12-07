#pragma once

#include <vector>
#include "DistributionFunction.hxx"

/*! @class CollisionOperator
 *
 *  @brief An abstract collision operator
 *
 *	This class defines the interface of a collision operator.
 *	It is designed to be inherited from and can not be instanciated.
 */
class CollisionOperator
{
public:
	CollisionOperator() = 0;
	virtual ~CollisionOperator();

	/*! @brief The operator computing the collision operation
	 *	
	 *	This operator defines the way the inherited classes will have to behave.
	 *
	 *	@param f The DistributionFunction to use
	 *
	 *	@param i The index on the x-axis of the point where the collision operation will be performed
	 *
	 *	@param j The index on the y-axis of the point where the collision operation will be performed
	 *
	 *	@returns A std::vector representing the new values at the given point of the distribution function after the collision step
	 */
	virtual std::vector<double> operator()(DistributionFunction const& f, size_t i, size_t j) const = 0;
};


/*! @class BGK
 *
 *  @brief A BGK collision operator
 *
 *	This class implements the abstract interface CollisionOperator.
 *	It represents the BGK model, which approximates Boltzmann's collision operator.
 */
class BGK : public CollisionOperator
{
protected:
	double m_nu;
	double m_gamma;
	double m_dt;

public:
	/*!	@brief The constructor of the class
	 *	
	 *	This constructor initializes the object and defines the fluid and simulation parameter that will be used by operator().
	 *
	 *	@param nu The kinematic viscosity of the fluid
	 *
	 *	@param gamma The adiabatic coefficient of the fluid
	 *
	 *	@param dt The time step of the simulation
	 */
	BGK(double nu, double gamma, double dt);

	/*! @brief The destructor of the class
	 *
	 */
	~BGK();

	/*! @brief The operator computing the collision operation
	 *	
	 *	This operator computes the collision step using the BGK approximation of Boltzmann's theoretical collision operator.
	 *	The computation uses a second order approximation when computing the equilibrium distribution function and
	 *	is therefore a good approximation of the Navier-Stokes model, as long as the Mach number stays low.
	 *	This operator uses the parameters defined at the construction of the instance.
	 *
	 *	@param f The DistributionFunction to use
	 *
	 *	@param i The index on the x-axis of the point where the collision operation will be performed
	 *
	 *	@param j The index on the y-axis of the point where the collision operation will be performed
	 *
	 *	@returns A std::vector representing the new values at the given point of the distribution function after the collision step
	 */
	virtual std::vector<double> operator()(DistributionFunction const& f, size_t i, size_t j) const;	
};