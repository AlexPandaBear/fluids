#pragma once

#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include "BoundaryCondition.hxx"

/*! @class LBMConfig
 *
 *	@todo Write the doc
 *	@todo set_boundary_conditions
 */
class LBMConfig
{
private:
	std::string m_sim_name;
	std::string m_sim_type; //RANS, URANS, LES, DNS ?
	std::string m_turb_model; //KE, KO, SA, None ?

	double m_mu;
	double m_lambda;
	double m_r;
	double m_gamma;

	size_t m_nb_dim; //1D, 2D, 3D ?
	std::vector<std::pair<double, double>> m_space_bounds;
	std::vector<size_t> m_space_grid;

	std::string m_velocity_scheme;

	std::pair<double, double> m_time_bounds; //Useless for RANS
	size_t m_nb_steps;

	std::vector<BoundaryCondition> m_bc;

	double m_output_dt; //Useless for RANS
	std::string m_output_type;

public:
	LBMConfig();
	LBMConfig(std::string config_file);
	~LBMConfig();

	void set_simulation_name(std::string name);
	void set_simulation_type(std::string type);
	void set_turbulence_model(std::string model);

	void set_fluid_properties(double mu, double lambda, double r, double gamma);

	void set_space(std::vector<std::pair<double, double>> bounds, std::vector<size_t> grid);
	void set_velocity_scheme(std::string scheme);
	void set_time(std::pair<double, double> bounds, size_t nb_steps);

	void set_output_parameters(double dt, std::string type);

	std::string get_simulation_name() const;
	std::string get_simulation_type() const;
	std::string get_turbulence_model() const;

	double get_mu() const;
	double get_lambda() const;
	double get_r() const;
	double get_gamma() const;

	size_t get_nb_dim() const;
	std::pair<double, double> get_space_bounds(size_t dimension) const;
	size_t get_grid(size_t dimension) const;
	std::vector<size_t> get_grid() const;

	std::string get_velocity_scheme() const;

	std::pair<double, double> get_time_bounds() const;
	size_t get_nb_steps() const;

	std::vector<BoundaryCondition> get_boundary_conditions() const;

	double get_output_dt() const;
	std::string get_output_type() const;

	void save(std::string file_name) const;
};