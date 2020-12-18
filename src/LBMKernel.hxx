#pragma once

#include <vector>
#include "DistributionFunction.hxx"
#include "Field.hxx"
#include "BoundaryCondition.hxx"
#include "CollisionOperator.hxx"
#include "StreamingOperator.hxx"

class LBMKernel
{
private:
	DistributionFunction m_f;

	std::vector<double> m_x, m_y;
	std::vector<double> m_time;

public:
	LBMKernel();
	~LBMKernel();

	void build_space(double x_min, double x_max, size_t nx, double y_min, double y_max, size_t ny);
	void build_time(double t_start, double t_end, size_t n);

	void simulate(CollisionOperator& collision, StreamingOperator& streaming, Field& f_out);
};