/**
 * @file osc_dym_more.cpp
 * @brief Verifies correctness of `libsolve_dym.a` (Harry Dym equation oscillatory
 *        solutions) via numerical differentiation, and demonstrates its usage.
 *
 * The program computes spatial/temporal derivatives (5‑point stencil) of the
 * library‑generated solutions and checks residuals against the PDE. It also
 * serves as a minimal example of how to invoke the library's functions.
 *
 * @author LI Ruomeng, liruomeng@zzu.edu.cn
 * @date 2026.02.24
 * @version 1.0
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <limits>
#include <cmath>
#include "solve_dym.h"

 // ============================================================================
 // Type Aliases and Using Declarations
 // ============================================================================

using std::abs;
using std::exp;
using std::log;
using std::max;
using std::vector;
using std::pair;
using std::array;

// ============================================================================
// Constants
// ============================================================================

namespace {
    constexpr int STENCIL_SIZE = 5;           ///< Size of the 5-point stencil
    constexpr int BOUNDARY_MARGIN = 2;        ///< Number of boundary points excluded on each side
}

// ----------------------------------------------------------------------------
// Stencil Coefficients for Finite Difference Approximations
// ----------------------------------------------------------------------------

/**
 * @brief Finite difference stencil coefficients for derivative approximations
 *
 * The 5-point stencil provides O(h⁴) accuracy for first and second derivatives:
 *
 * First derivative:  f'(x) ≈ [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)] / (12h)
 * Second derivative: f''(x) ≈ [-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)] / (12h²)
 *
 * Coefficient arrays:
 *   - d0: Identity stencil (for function values)
 *   - d1: First derivative stencil coefficients (numerator only)
 *   - d2: Second derivative stencil coefficients (numerator only)
 */
struct StencilCoefficients {
    real_t d0[STENCIL_SIZE] = { 0,  0,  1,  0,  0 };   ///< Identity: extracts center point
    real_t d1[STENCIL_SIZE] = { 1, -8,  0,  8, -1 };   ///< First derivative numerator
    real_t d2[STENCIL_SIZE] = { -1, 16,-30, 16, -1 };   ///< Second derivative numerator
};

// ----------------------------------------------------------------------------
// Core Functions
// ----------------------------------------------------------------------------

/**
 * @brief Computes spatial and temporal derivatives using 5-point stencil methods
 *
 * This function calculates the following derivatives for a function F(x,t) on a
 * uniform grid:
 *   - F00: Function values F(x,t)
 *   - F10: First spatial derivative ∂F/∂x
 *   - F20: Second spatial derivative ∂²F/∂x²
 *   - F01: First temporal derivative ∂F/∂t
 *
 * @tparam real_t Floating-point type (float, double, or long double)
 *
 * @param dx      Spatial step size (must be positive)
 * @param dt      Time step size (must be positive)
 * @param nx      Number of grid points in x-direction (must be ≥ 5)
 * @param nt      Number of grid points in t-direction (must be ≥ 5)
 * @param F       Input grid values of size nx × nt, column-major order:
 *                F[i + j*nx] = F(x_i, t_j) for i ∈ [0, nx-1], j ∈ [0, nt-1]
 * @param F00     Output array for function values on interior points
 * @param F10     Output array for ∂F/∂x on interior points
 * @param F20     Output array for ∂²F/∂x² on interior points
 * @param F01     Output array for ∂F/∂t on interior points
 *
 * @note All output arrays will be resized to (nx-4) × (nt-4) to accommodate
 *       the 5-point stencil which requires 2 boundary points on each side.
 *
 * @note Truncation error is O(Δx⁴) for spatial derivatives and O(Δt⁴) for
 *       temporal derivatives.
 *
 * @warning Undefined behavior if nx < 5 or nt < 5 (insufficient points for stencil)
 * @warning Undefined behavior if dx ≤ 0 or dt ≤ 0 (invalid step sizes)
 */
void compute_derivatives(
    real_t dx,
    real_t dt,
    size_t nx,
    size_t nt,
    const vector<real_t>& F,
    vector<real_t>& F00,
    vector<real_t>& F10,
    vector<real_t>& F20,
    vector<real_t>& F01)
{
    // -------------------------------------------------------------------------
    // Input Validation
    // -------------------------------------------------------------------------
    assert(dx > 0 && "Spatial step size dx must be positive");
    assert(dt > 0 && "Time step size dt must be positive");
    assert(nx >= STENCIL_SIZE && "Number of x-points must be at least 5");
    assert(nt >= STENCIL_SIZE && "Number of t-points must be at least 5");
    assert(F.size() == nx * nt && "Input array size mismatch");

    // -------------------------------------------------------------------------
    // Compute Interior Grid Dimensions
    // -------------------------------------------------------------------------
    // Exclude BOUNDARY_MARGIN points from each side
    const size_t nx_interior = nx - 2 * BOUNDARY_MARGIN;  // = nx - 4
    const size_t nt_interior = nt - 2 * BOUNDARY_MARGIN;  // = nt - 4
    const size_t interior_size = nx_interior * nt_interior;

    // Resize output arrays
    F00.resize(interior_size);
    F10.resize(interior_size);
    F20.resize(interior_size);
    F01.resize(interior_size);

    // -------------------------------------------------------------------------
    // Initialize Stencil Coefficients
    // -------------------------------------------------------------------------
    const StencilCoefficients stencil;

    // Precompute denominator factors for efficiency
    const real_t denom_dx = "12"_r * dx;
    const real_t denom_dx2 = "12"_r * dx * dx;
    const real_t denom_dt = "12"_r * dt;

    // -------------------------------------------------------------------------
    // Compute Derivatives at Interior Points
    // -------------------------------------------------------------------------
    for (int i = 0; i < nx_interior; ++i) {
        for (int j = 0; j < nt_interior; ++j) {
            real_t f00 = "0"_r;
            real_t f10 = "0"_r;
            real_t f20 = "0"_r;
            real_t f01 = "0"_r;

            // Apply 5-point stencil in both spatial and temporal directions
            for (size_t k = 0; k < STENCIL_SIZE; ++k) {
                for (size_t l = 0; l < STENCIL_SIZE; ++l) {
                    // Column-major indexing: index = i + j * nx
                    const size_t idx = i + k + (j + l) * nx;
                    const real_t f_val = F.at(idx);

                    // Accumulate weighted contributions for each derivative
                    f00 += stencil.d0[k] * stencil.d0[l] * f_val;
                    f10 += stencil.d1[k] * stencil.d0[l] * f_val;
                    f20 += stencil.d2[k] * stencil.d0[l] * f_val;
                    f01 += stencil.d0[k] * stencil.d1[l] * f_val;
                }
            }

            // Store results with appropriate scaling
            const size_t out_idx = i + j * nx_interior;
            F00.at(out_idx) = f00;
            F10.at(out_idx) = f10 / denom_dx;
            F20.at(out_idx) = f20 / denom_dx2;
            F01.at(out_idx) = f01 / denom_dt;
        }
    }
}

// ============================================================================
// Main Program - Verification Test
// ============================================================================

/**
 * @brief Main entry point for numerical verification
 *
 * This program verifies the accuracy of elliptic function solutions by:
 * 1. Computing solution values on a small grid around a test point
 * 2. Calculating numerical derivatives using 5-point stencils
 * 3. Evaluating residual errors against the governing PDEs
 * 4. Reporting the maximum absolute error (epsilon)
 *
 * The governing equations being verified are:
 *   - X_x = S^2  (first residual: epsilon_1)
 *   - X_t = S S_{xx} - 2 S_{x}^2 + (2-m) S^2
 * With R = S^2 as a function of (X, T) rather than (x, t), R(X, T) satisfies
 * the Harry Dym equation:
 *   R_T = (1/2) * R^3 * R_{XXX}
 *
 * @param argc Argument count (unused)
 * @param args Argument values (unused)
 * @return 0 on success, non-zero on failure
 */
int main(int argc, char* args[])
{
    (void)argc;  // Suppress unused parameter warning
    (void)args;  // Suppress unused parameter warning

    // -------------------------------------------------------------------------
    // Solution Parameters
    // -------------------------------------------------------------------------

    /// Number of waves (cuspons + solitons) in the solution (positive integer)
    constexpr size_t N = 3;

    /// Number of cuspon components (0 <= Np <= N), 
    /// and number of soliton components (N - Np)
    constexpr size_t Np = 0;

    /// Elliptic modulus parameter (0 < m < 1)
    const real_t m = "0.7"_r;

    // -------------------------------------------------------------------------
    // Grid Parameters (you CAN choose your own parameter)
    // -------------------------------------------------------------------------

    /// Spatial step size
    const real_t dx = "0.00369"_r;

    /// Temporal step size
    const real_t dt = "0.0056789"_r;

    /// Center point in space
    const real_t x0 = "0.3456789"_r;

    /// Center point in time
    const real_t t0 = "0.38324"_r;

    const int nx = 2004;
    const int nt = 2004;
    // Grid points for 5-point stencil (2 points on each side of center)
    vector<real_t> x(nx);
    for (size_t i = 0; i < nx; ++i)
    {
        x[i] = x0 + i * dx - 2 * dx;
    }

    vector<real_t> t(nt);
    for (size_t i = 0; i < nt; ++i)
    {
        t[i] = t0 + i * dt - 2 * dt;
    }

    // -------------------------------------------------------------------------
    // Solution Parameters (Spectral Data) (you CAN choose your own parameter)
    // -------------------------------------------------------------------------

    /// Spectral parameters z_j
    /// First Np values and last (N-Np) values should be distinct
    const real_t dz = "0.01"_r;
    const real_t z0 = "0.2"_r; 
    const real_t z1 = "0.3"_r;
    vector<real_t> z(N);
    for (size_t i=0; i < Np; ++i)
    {
        z[i] = z0 + i * dz;
    }
    for (size_t i=0; i < N - Np; ++i)
    {
        z[Np + i] = z1 + i * dz;
    }

    /// Phase parameters
    const vector<real_t> omega (N, "0"_r);

    // -------------------------------------------------------------------------
    // Compute Solution on Grid (Modify cautiously)
    // -------------------------------------------------------------------------

    /// Solution arrays: S, X, T
    vector<real_t> S, X, T;

    compute_solution(x, t, N, Np, z, omega, m, S, X, T);

    // -------------------------------------------------------------------------
    // Compute Numerical Derivatives
    // -------------------------------------------------------------------------

    /// Derivative arrays for S (subscript: 00=value, 10=∂x, 20=∂²x, 01=∂t)
    vector<real_t> S00, S10, S20, S01;
    vector<real_t> X00, X10, X20, X01;
    vector<real_t> T00, T10, T20, T01;

    compute_derivatives(dx, dt, nx, nt, S, S00, S10, S20, S01);
    compute_derivatives(dx, dt, nx, nt, X, X00, X10, X20, X01);
    compute_derivatives(dx, dt, nx, nt, T, T00, T10, T20, T01);

    // -------------------------------------------------------------------------
    // Compute Residual Errors
    //-------------------------------------------------------------------------

    /// Maximum absolute error across all interior points
    real_t epsilon = -std::numeric_limits<real_t>::infinity();

    for (size_t i = 0; i < X00.size(); ++i) {
        // Residual 1: X_x - S^2 = 0
        const real_t epsilon_1 = abs(X10.at(i) - S00.at(i) * S00.at(i));

        // Residual 2: X_t - S S_{xx} + 2 S_x^2 - (2-m) S^2 = 0
        const real_t epsilon_2 = abs(
            X01.at(i)
            - S00.at(i) * S20.at(i)
            + "2"_r * S10.at(i) * S10.at(i)
            - ("2"_r - m) * S00.at(i) * S00.at(i)
        );

        // Track maximum error
        epsilon = max(epsilon, epsilon_1);
        epsilon = max(epsilon, epsilon_2);
    }

    // -------------------------------------------------------------------------
    // Output Results
    // -------------------------------------------------------------------------
    std::cout << std::scientific
        << std::setprecision(6)
        << "epsilon = "
        << std::setw(18)
        << epsilon
        << std::endl;

    return (epsilon < 1e-6) ? 0 : 1;  // Return success if error is below tolerance
}