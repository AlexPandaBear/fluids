#pragma once

#include <vector>
#include "DistributionFunction.hxx"

/*! @class CollisionOperator
 *
 *  @brief An abstract collision operator
 *
 *	This class defines the interface of a collision operator.
 *	It is designed to be inherited from and can not be instanciated.
 *
 *	@todo Change df API
 */
class CollisionOperator
{
public:
	//CollisionOperator() = 0;
	virtual ~CollisionOperator();

	/*! @brief The operator computing the collision operation
	 *	
	 *	This operator defines the way the inherited classes will have to behave.
	 *
	 *	@param f The DistributionFunction to use
	 *
	 *	@param dt The length of the time step
	 */
	virtual void operator()(DistributionFunction& f, double dt) const;
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
	double m_mu;
	double m_gamma;

public:
	/*!	@brief The constructor of the class
	 *	
	 *	This constructor initializes the object and defines the fluid and simulation parameter that will be used by operator().
	 *
	 *	@param mu The viscosity of the fluid
	 *
	 *	@param gamma The adiabatic coefficient of the fluid
	 */
	BGK(double mu, double gamma);

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
	 *	@param dt The length of the time step
	 */
	virtual void operator()(DistributionFunction& f, double dt) const;	
};