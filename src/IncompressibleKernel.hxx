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
	const double m_theta;

	const double m_Lx, m_Ly;
	const size_t m_nx, m_ny;
	const size_t m_nb_pts;
	const double m_dx, m_dy;
	const double m_dxdy, m_dtdx, m_dtdy;

	const double m_rho;
	const double m_mu;
	const double m_nu;

	const double m_accuracy_u, m_accuracy_v, m_accuracy_p;
	const bool m_clean_pressure;

	std::vector<std::vector<bool>> m_body;

	const double m_nu_x, m_nu_y;
	const double m_delta_x, m_delta_y;

	std::vector<std::vector<double>> m_U, m_V, m_P;
	std::vector<std::vector<double>> U_tilde, V_tilde;

	std::vector<double> diag_alpha, diag_beta1, diag_beta2, diag_gamma1, diag_gamma2;
	std::vector<double> diag_delta1u, diag_delta2u,	diag_delta1v, diag_delta2v;
	std::vector<double> diag_Au, diag_B1u, diag_B2u, diag_D1, diag_D2, diag_D3, diag_D4;
	std::vector<double> diag_Av, diag_C1v, diag_C2v, diag_E1, diag_E2, diag_E3, diag_E4;
	std::vector<double> diag_A, diag_B1, diag_B2, diag_C1, diag_C2;

	std::vector<double> vec_X, vec_B;
/*
	void getPressureFieldAt(size_t t);
	void getVelocityFieldsAt(size_t t);
*/
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

//	void updateVelocityField(size_t t);
	void computePressureField(PoissonSolver& solver, size_t t);

	void printSimulationProgression(size_t t) const;

	void loadFieldsAt(size_t t);
	void saveFieldsAt(size_t t) const;

	void computeMatrixA();
	void computeVectorB();

	double safe_access_VectorX(size_t block, size_t i, size_t j) const;

	double GS_iteration_U();
	double GS_iteration_V();
	double GS_iteration_P();
	void solve_GaussSeidel();

public:
	IncompressibleKernel(DataKeeper& data, double tMax, size_t nb_steps, double theta, double Lx, double Ly, size_t nx, size_t ny, double rho, double mu, double accuracy_u, double accuracy_v, double accuracy_p, bool clean_pressure);
	~IncompressibleKernel();

	void defineBoundaryConditions(std::vector<double> v_U_bc, std::vector<double> v_V_bc, std::vector<double> v_P_bc);
	void defineBody(std::vector<std::vector<bool>> v_body);

//	void simulate();
	void simulate_theta();
};