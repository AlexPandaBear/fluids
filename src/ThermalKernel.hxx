#pragma once

#include <vector>
#include <iostream>

#include "DataKeeper.hxx"

class ThermalKernel
{
private:
	DataKeeper& m_data;

	double m_tMax;
	size_t m_nb_steps;
	double m_dt;
	double m_theta;
	double m_accuracy;

	double m_Lx, m_Ly;
	size_t m_nx, m_ny;
	double m_dx, m_dy;
	size_t m_nb_eq;

	double m_lambda;
	double m_rho;
	double m_cp;
	double m_K;

	double m_T_BC;

	double m_kx, m_ky;
	std::vector<std::vector<double>> m_u, m_v, m_T;

	std::vector<double> mat_A_a, mat_A_b, mat_A_c, mat_A_d, mat_A_e;
	std::vector<double> vec_X, vec_B;

	void display();

	void computeMatrixA(size_t t);
	void computeVectorB(size_t t);
	void solve_GaussSeidel();

	void printSimulationProgression(size_t t) const;

public:
	ThermalKernel(DataKeeper& data, double tMax, size_t nb_steps, double theta, double accuracy, double Lx, double Ly, size_t nx, size_t ny, double lambda, double rho, double cp, double T_BC);
	~ThermalKernel();
	
	void simulate();
};