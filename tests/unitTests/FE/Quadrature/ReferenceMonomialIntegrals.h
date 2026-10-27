// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVMP_TEST_REFERENCE_MONOMIAL_INTEGRALS_H
#define SVMP_TEST_REFERENCE_MONOMIAL_INTEGRALS_H

/** @file ReferenceMonomialIntegrals.h
 * @brief Independent analytic moments for quadrature tests only.
 */
namespace svmp::FE::quadrature::test {

/** @brief Integrate x^power on [-1,1].
 * @throws std::invalid_argument If power is negative.
 */
long double line_monomial_integral(int power);

/** @brief Integrate x^a y^b on the unit right triangle using bounded ratios.
 * @throws std::invalid_argument If either exponent is negative.
 */
long double triangle_monomial_integral(int a, int b);

} // namespace svmp::FE::quadrature::test

#endif // SVMP_TEST_REFERENCE_MONOMIAL_INTEGRALS_H
