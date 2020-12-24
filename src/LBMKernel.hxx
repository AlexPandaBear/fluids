#pragma once

#include <vector>
#include <string>
#include <utility>
#include "LBMConfig.hxx"
#include "DistributionFunction.hxx"
#include "Field.hxx"
#include "BoundaryCondition.hxx"
#include "CollisionOperator.hxx"
#include "StreamingOperator.hxx"

/*! @class LBMKernel
 *
 *	@todo Write the doc
 */
class LBMKernel
{
private:
	LBMConfig m_config;

	DistributionFunction m_f;

	CollisionOperator m_collision;
	StreamingOperator m_streaming;

	std::vector<std::vector<double>> m_space;
	std::vector<double> m_time;

	void build_space();
	void build_time();
	void initialize_f();
	void initialize_operators();

public:
	LBMKernel(LBMConfig& config);
	LBMKernel(std::string config_file);
	~LBMKernel();

	/*! @todo Finish API change
	 */
	void set_initial_state(Field& rho, Field& U, Field& E);

	void simulate(Field& f_out);
};