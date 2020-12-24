#pragma once

#include <vector>
#include "DistributionFunction.hxx"
#include "BoundaryCondition.hxx"
#include "Field.hxx"

/*! @class StreamingOperator
 *
 *	This class defines the interface of a streaming operator.
 *	It is designed to be inherited from and can not be instanciated.
 *
 *	@todo Change df API
 */
class StreamingOperator
{
public:
	//StreamingOperator(std::vector<BoundaryCondition>& BCs) = 0;
	virtual ~StreamingOperator();

	/*! @brief The method making sure the streaming operator is ready for computation
	 *
	 * 	@param f The DistributionFunction that will be used for computation
	 */
	virtual void initialize(DistributionFunction const& f);

	/*! @brief The operator computing the streaming operation
	 *	
	 *	This operator defines the way the inherited classes will have to behave.
	 *
	 *	@param f The DistributionFunction to stream
	 */
	virtual void operator()(DistributionFunction& f);
};

/*! @class Stream2D
 *
 *	This class implements the abstract interface StreamingOperator.
 *	As its name hints, it was designed for 2D simulations.
 */
class Stream2D : public StreamingOperator
{
private:
	size_t m_nx, m_ny, m_nv;
	Field m_tmp;

	size_t m_nb_BC;
	std::vector<BoundaryCondition> m_BC;

public:
	/*! @brief The constructor of the class 
	 *	
	 *	This constructor prepares the instance to perform the streaming operation on a specific simulation.
	 *	It initializes all its member attributes to the right values given the shape of the DistributionFunction that will be handled.
	 *
	 *	@param BC The boundary conditions that will be used
	 */
	Stream2D(std::vector<BoundaryCondition> BC);

	/*! @brief The destructor of the class
	 */
	virtual ~Stream2D();

	/*! @brief The method making sure the streaming operator is ready for computation
	 *
	 * 	@param f The DistributionFunction that will be used for computation
	 */
	virtual void initialize(DistributionFunction const& f);

	/*! @brief The operator computing the streaming operation
	 *
	 *	This operator performs the streaming step on the given DistributionFunction object.
	 *	
	 *	@param f The distribution function to stream
	 *
	 *	@warning Using this operator will modify f
	 */
	virtual void operator()(DistributionFunction& f);
};