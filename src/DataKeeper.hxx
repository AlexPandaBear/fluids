#pragma once

#include <memory>
#include <string>
#include <fstream>
#include <iostream>

/**
 * A class responsible for the storage of the simulation data and the management of the simulation files
 */
class DataKeeper
{
private:
	size_t m_nb_steps;
	size_t m_nx;
	size_t m_ny;

	double m_coef_1;
	double m_coef_2;

	std::unique_ptr<double[]> ptr_T;
	std::unique_ptr<double[]> ptr_P;
	std::unique_ptr<double[]> ptr_U;
	std::unique_ptr<double[]> ptr_V;

public:
	/**
	 * The constructor of the class
	 */
	DataKeeper();

	/**
	 * The destructor of the class
	 */
	~DataKeeper();

	/**
	 * A method reshaping the containers to fit the data computed by a new simulation
	 * 
	 * @param nb_steps The number of steps in the temporal integration
	 *
	 * @param nx The number of points of computation along the x-axis
	 *
	 * @param ny The number of points of computation along the y-axis
	 *
	 * @warning When calling this method, any data previously stored in this instance will be lost
	 */
	void reset_size(size_t nb_steps, size_t nx, size_t ny);

	/**
	 * A method accessing the temperature value stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 */
	double getTemperatureAt(size_t t, size_t i, size_t j) const;

	/**
	 * A method accessing the pressure value stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 */
	double getPressureAt(size_t t, size_t i, size_t j) const;

	/**
	 * A method accessing the value of the x-coordinate of the velocity stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 */
	double getXVelocityAt(size_t t, size_t i, size_t j) const;

	/**
	 * A method accessing the value of the y-coordinate of the velocity stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 */
	double getYVelocityAt(size_t t, size_t i, size_t j) const;

	/**
	 * A method setting the temperature value stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 *
	 * @param T The value to store
	 *
	 * @warning Calling this method will overwrite the previously stored value
	 */
	void setTemperatureAt(size_t t, size_t i, size_t j, double T);

	/**
	 * A method setting the pressure value stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 *
	 * @param P The value to store
	 *
	 * @warning Calling this method will overwrite the previously stored value
	 */
	void setPressureAt(size_t t, size_t i, size_t j, double P);

	/**
	 * A method setting the value of the x-coordinate of the velocity stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 *
	 * @param U The value to store
	 *
	 * @warning Calling this method will overwrite the previously stored value
	 */
	void setXVelocityAt(size_t t, size_t i, size_t j, double U);

	/**
	 * A method setting the value of the y-coordinate of the velocity stored at a specific moment and location
	 *
	 * @param t The time index of the value
	 *
	 * @param i The x-coordinate index of the point
	 *
	 * @param j The y-coordinate index of the point
	 *
	 * @param V The value to store
	 *
	 * @warning Calling this method will overwrite the previously stored value
	 */
	void setYVelocityAt(size_t t, size_t i, size_t j, double V);

	/**
	 * A method saving all the current state in the instance in a file in order to get back to this state later
	 *
	 * @param file_name The full name (with the eventual path) of the file to create
	 *
	 * @warning If the name chosen matches an existing file, this method will overwrite the existing file
	 */
	void saveData(std::string file_name) const;

	/**
	 * A method loading all the data from a file in order to resume a previous state of an instance of the class
	 *
	 * @param file_name The full name (with the eventual path) of the file to load
	 */
	void loadData(std::string file_name);
};