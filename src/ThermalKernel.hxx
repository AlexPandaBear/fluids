#pragma once

#include <vector>
#include <iostream>

#include "DataKeeper.hxx"

/**
 * A class designed to be called by an instance of SimManager in order to solve the heat equation.
 */

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
	/**
	 * The constructor of the class
	 *
	 * @param data The DataKeeper object in which the computed solution will be stored
	 *
	 * @param tMax The duration of the simulation
	 *
	 * @param nb_steps The number of steps to perform for the temporal integration
	 *
	 * @param theta The parameter of the theta-scheme used for the temporal integration
	 *
	 * @param accuracy The accuracy demanded for the computed temperatures (stop condition for the Gauss-Seidel algorithm)
	 *
	 * @param Lx The length along the x-axis of the domain to simulate
	 *
	 * @param Ly The length along the y-axis of the domain to simulate
	 *
	 * @param nx The number of points of computation along the x-axis
	 *
	 * @param ny The number of points of computation along the y-axis
	 *
	 * @param lambda The thermal conductivity of the fluid
	 *
	 * @param rho The volumic mass of the fluid
	 *
	 * @param cp The heat capacity per unit of mass of the fluid
	 *
	 * @param T_BC The temperature at the boundary of the domain
	 */
	ThermalKernel(DataKeeper& data, double tMax, size_t nb_steps, double theta, double accuracy, double Lx, double Ly, size_t nx, size_t ny, double lambda, double rho, double cp, double T_BC);
	
	/**
	 * The destructor of the class
	 */
	~ThermalKernel();
	
	/**
	 * A method solving the heat equation thanks to the parameters previously defined
	 *
	 * @warning To work properly, this method has to be called after the computation of the flow.
	 *
	 * @warning To work properly, this method has to be called after defining the initial temperature field in the DataKeeper object used to create the instance.
	 */
	void simulate();
};