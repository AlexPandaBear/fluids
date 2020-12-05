#pragma once

#include <functional>
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
protected:
	DistributionFunction const& m_f;

public:
	CollisionOperator() = 0;
	virtual ~CollisionOperator();

	/*! @brief A setter for the distribution function
	 *  
	 *  This method sets the member DistributionFunction object which will be used for computation.
	 *	This object is stored by reference to avoid useless copies.
	 *
	 *  @param f The DistributionFunction object to link to this instance
	 */
	void set_f(DistributionFunction const& f);
	
	virtual double operator()(size_t i, size_t j, size_t k) const = 0;
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
	double m_r;

	field..

	std::function<double(size_t, size_t, size_t, field, field, field, field ..)> m_f_eq;

public:
	BGK(double nu, double r, std::function<> f_eq);
	~BGK();
	
	virtual double operator()(size_t i, size_t j, size_t k) const;	
};