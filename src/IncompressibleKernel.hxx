#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include "DataKeeper.hxx"

class IncompressibleKernel
{
private:
	DataKeeper &m_data;

	double m_tMax;
	size_t m_nb_steps;
	double m_dt;

	double m_Lx, m_Ly;
	size_t m_nx, m_ny;
	double m_dx, m_dy;

	std::vector<std::vector<bool>> v_body;

	double m_rho;
	double m_mu;
	double m_nu;

	double m_err_max;

	std::vector<std::vector<double>> m_U;
	std::vector<std::vector<double>> m_V;
	std::vector<std::vector<double>> m_P;
	std::vector<std::vector<double>> m_P_tmp;

	double m_coef_1;
	double m_coef_2;
	double m_coef_3;

	void initialize();
	void getAllFieldsAt(size_t t);
	void getVelocityFieldsAt(size_t t);

	double dudx(size_t i, size_t j) const;
	double dudy(size_t i, size_t j) const;
	double dvdx(size_t i, size_t j) const;
	double dvdy(size_t i, size_t j) const;

	double laplacian_u(size_t i, size_t j) const;
	double laplacian_v(size_t i, size_t j) const;

	double advection_u(size_t i, size_t j) const;
	double advection_v(size_t i, size_t j) const;
	double diffusion_u(size_t i, size_t j) const;
	double diffusion_v(size_t i, size_t j) const;
	double pressure_u(size_t i, size_t j) const;
	double pressure_v(size_t i, size_t j) const;

	void updateVelocityField(size_t t);
	void updatePressureField(size_t t, double mean_error_max);

	void printSimulationProgression(size_t t) const;

public:
	IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu, double err_max);
	~IncompressibleKernel();

	void simulate();
};