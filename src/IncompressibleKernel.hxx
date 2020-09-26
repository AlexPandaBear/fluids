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

	const double m_tMax;
	const size_t m_nb_steps;
	const double m_dt;

	const double m_Lx, m_Ly;
	const size_t m_nx, m_ny;
	const double m_dx, m_dy;

	const double m_rho;
	const double m_mu;
	const double m_nu;

	std::vector<std::vector<double>> m_U;
	std::vector<std::vector<double>> m_V;
	std::vector<std::vector<double>> m_P;

	void getPressureFieldAt(size_t t);
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
	void computePressureField(PoissonSolver& solver, size_t t);

	void printSimulationProgression(size_t t) const;

public:
	IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu);
	~IncompressibleKernel();

	void defineBoundaryConditions(std::vector<double> v_U_bc, std::vector<double> v_V_bc, std::vector<double> v_P_bc);

	void simulate();
};