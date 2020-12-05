#pragma once

#include <vector>

/*! @brief A LBM velocity scheme
 *	
 *	This class represents a LBM DxQy velocity scheme, x being the number dimensions of the fluid domain and y the number of discrete velocities.
 */
class VelocityScheme
{
private:
	size_t m_nb_dim, m_nb_vel;

	std::vector<double*> m_velocities;
	double* m_weights;

public:
	VelocityScheme(size_t nb_dim, size_t nb_vel);
	~VelocityScheme();
	
};