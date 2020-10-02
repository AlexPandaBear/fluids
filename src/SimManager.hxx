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
	double m_cp;
	double m_mu;

	double m_theta_inc;
	double m_accuracy_u, m_accuracy_v, m_accuracy_p;
	bool m_clean_pressure;

	double m_theta_th;
	double m_accuracy_th;

	std::vector<double> m_U_BC, m_V_BC, m_P_BC;

	std::vector<std::vector<bool>> m_body;

	bool m_temp_BC, m_flux_BC;
	double m_T_BC;

	double m_press_err_max = 10.;

	DataKeeper m_data;

public:
	SimManager();
	~SimManager();

	void defineTimeParameters(double tMax, size_t nb_steps);
	void defineGridParameters(double Lx, double Ly, size_t nx, size_t ny);
	void defineFluidProperties(double lambda, double rho, double cp, double mu);
	void defineBody(std::vector<std::vector<bool>> body);
	void defineInitialState(std::vector<std::vector<double>> U0, std::vector<std::vector<double>> V0, std::vector<std::vector<double>> T0);
	void defineUniformDynamicBoundaryConditions(double U_bc, double V_bc, double P_bc);
	void defineDynamicBoundaryConditions(std::vector<double> v_U_bc, std::vector<double> v_V_bc, std::vector<double> v_P_bc);
	void defineThermalBoundaryConditions(std::string type, double value);
	void defineFlowIntegrationParameters(double theta, double accuracy_u, double accuracy_v, double accuracy_p, bool clean_pressure);
	void defineThermalIntegrationParameters(double theta, double accuracy);

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
};