#pragma once

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>
#include "IncompressibleKernel.hxx"

class SimManager
{
private:
	double m_tMax;
	size_t m_nb_steps;
	double m_dt;

	double m_Lx, m_Ly;
	size_t m_nx, m_ny;
	double m_dx, m_dy;

	double m_lambda;
	double m_rho;
	double m_cv;
	double m_mu;
	double m_nu;
	double m_Cv;
	double m_coef_x, m_coef_y;

	bool m_temp_BC, m_flux_BC;
	double m_BC_value;

	std::vector<std::vector<std::vector<double>>> m_T;
	std::vector<std::vector<std::vector<double>>> m_U;
	std::vector<std::vector<std::vector<double>>> m_V;
	std::vector<std::vector<std::vector<double>>> m_P;

	double diffusion(size_t t, size_t i, size_t j) const;
	double advection(size_t t, size_t i, size_t j) const;
	double creation(size_t t, size_t i, size_t j) const;

	void computeFlow();

public:
	SimManager();
	~SimManager();

	void defineTimeParameters(double tMax, size_t nb_steps);
	void defineGridParameters(double Lx, double Ly, size_t nx, size_t ny);
	void defineFluidProperties(double lambda, double rho, double cv, double mu);
	void defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> P0, std::vector<std::vector<double>> T0);
	void defineBoundaryConditions(std::string type, double value);
	void launchSimulation();

	std::vector<std::vector<std::vector<double>>> getResults() const;
};