#pragma once

#include <vector>
#include <stdexcept>
#include <utility>
#include "Field.hxx"

class BoundaryCondition
{
public:
	//BoundaryCondition() = 0;
	virtual ~BoundaryCondition();
	
	virtual void apply(Field& f) const;
};

class D2Q9BounceBack : public BoundaryCondition
{
protected:
	size_t m_nb_nodes;
	std::vector<std::pair<size_t, size_t>> m_ij;
	std::vector<std::vector<size_t>> m_dir;

public:
	D2Q9BounceBack(std::vector<std::pair<size_t, size_t>> const& ij, std::vector<std::vector<size_t>> const& dir);
	virtual ~D2Q9BounceBack();

	virtual void apply(Field& f) const;
};