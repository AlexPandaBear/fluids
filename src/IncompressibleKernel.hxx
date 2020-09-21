#pragma once

#include <vector>
#include <cmath>

class IncompressibleKernel
{
private:
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

	std::vector<std::vector<std::vector<double>>> m_U;
	std::vector<std::vector<std::vector<double>>> m_V;
	std::vector<std::vector<std::vector<double>>> m_P;

	void initialize();

	double dudx(size_t t, size_t i, size_t j) const;
	double dudy(size_t t, size_t i, size_t j) const;
	double dvdx(size_t t, size_t i, size_t j) const;
	double dvdy(size_t t, size_t i, size_t j) const;

	double laplacian_u(size_t t, size_t i, size_t j) const;
	double laplacian_v(size_t t, size_t i, size_t j) const;

	double advection_u(size_t t, size_t i, size_t j) const;
	double advection_v(size_t t, size_t i, size_t j) const;
	double diffusion_u(size_t t, size_t i, size_t j) const;
	double diffusion_v(size_t t, size_t i, size_t j) const;
	double pressure_u(size_t t, size_t i, size_t j) const;
	double pressure_v(size_t t, size_t i, size_t j) const;

	void updateVelocityField(size_t t);
	void updatePressureField(size_t t, double mean_error_max);

public:
	IncompressibleKernel();
	~IncompressibleKernel();

	void defineTimeParameters(double tMax, size_t nb_steps);
	void defineGridParameters(double Lx, double Ly, size_t nx, size_t ny);
	void defineBodyShape(std::vector<std::vector<bool>> body_shape);
	void defineFluidProperties(double rho, double mu);
	void defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> P0);
	void simulate();
	std::vector<std::vector<std::vector<double>>> getU() const;
	std::vector<std::vector<std::vector<double>>> getV() const;
	std::vector<std::vector<std::vector<double>>> getP() const;
};