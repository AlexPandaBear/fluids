#pragma once

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>
#include "IncompressibleKernel.hxx"
#include "ThermalKernel.hxx"
#include "DataKeeper.hxx"

class SimManager
{
private:
	double m_tMax;
	size_t m_nb_steps;

	double m_Lx, m_Ly;
	size_t m_nx, m_ny;

	double m_lambda;
	double m_rho;
	double m_cv;
	double m_mu;

	bool m_temp_BC, m_flux_BC;
	double m_BC_value;

	double m_press_err_max = 10.;

	DataKeeper m_data;

public:
	SimManager();
	~SimManager();

	void defineTimeParameters(double tMax, size_t nb_steps);
	void defineGridParameters(double Lx, double Ly, size_t nx, size_t ny);
	void defineFluidProperties(double lambda, double rho, double cv, double mu);
	void defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> P0, std::vector<std::vector<double>> T0);
	void defineBoundaryConditions(std::string type, double value);
	void defineErrorTolerance(double press_err_max);

	void launchSimulation();

	double getTemperatureAt(size_t t, size_t i, size_t j) const;
	double getPressureAt(size_t t, size_t i, size_t j) const;
	double getXVelocityAt(size_t t, size_t i, size_t j) const;
	double getYVelocityAt(size_t t, size_t i, size_t j) const;

	void setTemperatureAt(size_t t, size_t i, size_t j, double T);
	void setPressureAt(size_t t, size_t i, size_t j, double P);
	void setXVelocityAt(size_t t, size_t i, size_t j, double U);
	void setYVelocityAt(size_t t, size_t i, size_t j, double V);

	void saveData(std::string file_name) const;
	void loadData(std::string file_name);
/*
	double get_tMax() const;
	size_t get_nb_steps() const;
	double get_dt() const;
	double get_Lx() const;
	double get_Ly() const;
	size_t get_nx() const;
	size_t get_ny() const;
	double get_dx() const;
	double get_dy() const;
	double get_lambda() const;
	double get_rho() const;
	double get_cv() const;
	double get_mu() const;
	double get_nu() const;
	double get_Cv() const;
	bool get_temp_BC() const;
	bool get_flux_BC() const;
	double get_BC_value() const;
	double get_press_err_max() const;
*/
};