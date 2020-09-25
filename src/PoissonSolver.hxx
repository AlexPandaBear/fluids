#pragma once

#include <vector>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <iostream>

class PoissonSolver
{
private:
	const double m_Lx, m_Ly;
	const size_t m_nx, m_ny;
	const double m_dx, m_dy;

	const size_t m_nb_eq;

	const double m_coef_1;
	const double m_coef_2;

	std::vector<std::vector<double>> mat_lower;
	std::vector<std::vector<double>> mat_upper;

	std::vector<double> v_x;
	std::vector<double> v_y;
	std::vector<double> v_f;

	double matrix_coef(size_t i, size_t k);
	void generate_LU();
	void solve_for_y();
	void solve_for_x();

public:
	PoissonSolver(double Lx, double Ly, size_t nx, size_t ny);
	~PoissonSolver();

	void definePoissonFunction(std::vector<double> v_values);
	void defineBoundaryConditions(double uniform_value);
	void defineBoundaryConditions(std::vector<double> v_up, std::vector<double> v_down, std::vector<double> v_left, std::vector<double> v_right);
	std::vector<double> solve();
};