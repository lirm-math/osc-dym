#pragma once
#include <vector>
#include <limits>
#include "real_t.h"

void compute_solution(
	const std::vector<real_t>& x,
	const std::vector<real_t>& t,
	size_t N, size_t Np,
	const std::vector<real_t>& z,
	const std::vector<real_t>& omega,
	real_t m,
	std::vector<real_t>& S,
	std::vector<real_t>& X,
	std::vector<real_t>& T);
