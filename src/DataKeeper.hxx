#pragma once

#include <memory>
#include <string>
#include <fstream>

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
	DataKeeper();
	~DataKeeper();

	void reset_size(size_t nb_steps, size_t nx, size_t ny);

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