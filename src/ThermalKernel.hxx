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

	double m_Lx, m_Ly;
	size_t m_nx, m_ny;
	double m_dx, m_dy;

	double m_lambda;
	double m_rho;
	double m_cv;
	double m_Cv;
	double m_coef_x, m_coef_y;

	std::vector<std::vector<double>> m_T;
	std::vector<std::vector<double>> m_U;
	std::vector<std::vector<double>> m_V;

	bool m_temp_BC, m_flux_BC;
	double m_BC_value;

	void getAllFieldsAt(size_t t);

	double diffusion(size_t i, size_t j) const;
	double advection(size_t i, size_t j) const;
	double creation(size_t i, size_t j) const;

	void printSimulationProgression(size_t t) const;

public:
	ThermalKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double lambda, double rho, double cv, bool temp_BC, bool flux_BC, double BC_value);
	~ThermalKernel();
	
	void simulate();
};