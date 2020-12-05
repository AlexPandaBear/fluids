#pragma once

#include <vector>
#include <cmath>
#include "DataKeeper.hxx"

class DataProcessor
{
private:
	bool m_ready;

	DataKeeper &m_data;
	
	size_t m_nx, m_ny, m_nb_steps;
	double m_dx, m_dy, m_dt;

	double dudx(size_t t, size_t i, size_t j) const;
	double dvdy(size_t t, size_t i, size_t j) const;

public:
	DataProcessor(DataKeeper& data);
	~DataProcessor();

	bool is_ready() const;
	void define_parameters(DataKeeper &data, size_t nx, size_t ny, size_t nb_steps, double dx, double dy, double dt);
	
	std::vector<std::vector<double>> norm_V_field(size_t t) const;
	std::vector<std::vector<double>> div_V_field(size_t t) const;
};