#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include "DataKeeper.hxx"
#include "PoissonSolver.hxx"

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

	double m_rho;
	double m_mu;
	double m_nu;

	std::vector<std::vector<double>> m_U;
	std::vector<std::vector<double>> m_V;
	std::vector<std::vector<double>> m_P;

	void getAllFieldsAt(size_t t);
	void getVelocityFieldsAt(size_t t);

	double dudx(size_t i, size_t j) const;
	double dudy(size_t i, size_t j) const;
	double dvdx(size_t i, size_t j) const;
	double dvdy(size_t i, size_t j) const;

	double d2udx2(size_t i, size_t j) const;
	double d2udy2(size_t i, size_t j) const;
	double d2vdx2(size_t i, size_t j) const;
	double d2vdy2(size_t i, size_t j) const;

	double d2uvdxdy(size_t i, size_t j) const;	
	double d2u2dx2(size_t i, size_t j) const;
	double d2v2dy2(size_t i, size_t j) const;

	double laplacian_u(size_t i, size_t j) const;
	double laplacian_v(size_t i, size_t j) const;

	double advection_u(size_t i, size_t j) const;
	double advection_v(size_t i, size_t j) const;
	double diffusion_u(size_t i, size_t j) const;
	double diffusion_v(size_t i, size_t j) const;
	double pressure_u(size_t i, size_t j) const;
	double pressure_v(size_t i, size_t j) const;

	void updateVelocityField(size_t t);
	void updatePressureField(PoissonSolver& solver, size_t t);

	void printSimulationProgression(size_t t) const;

public:
	IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu, double err_max);
	~IncompressibleKernel();

	void simulate();
};