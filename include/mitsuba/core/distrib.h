#pragma once
#if !defined(__DISTRIB_H)
#define __DISTRIB_H

#include <mitsuba/mitsuba.h>
#include <vector>

MTS_NAMESPACE_BEGIN

/**
 * @brief Implements a 1D discrete probability distribution for importance sampling.
 *
 * This class takes a 1D array of function values (e.g., from a piece-wise
 * constant function) and builds a data structure that allows for efficient
 * sampling. Samples are drawn proportionally to the function's values.
 *
 * The implementation is based on the inversion method. It precomputes the
 * cumulative distribution function (CDF) to enable fast sampling.
 */
struct MTS_EXPORT_CORE DistribTable {
	DistribTable() = default;

	/**
	 * @brief Constructs a 1D distribution from a given array of weights.
	 * @param f A pointer to an array of 'n' floating-point values.
	 * @param n The number of values in the array 'f'.
	 */
	DistribTable(const Float *f, uint32_t n);

	/**
	 * @brief Draws a sample from the discrete distribution.
	 * @param u A uniformly distributed random number in [0, 1).
	 * @param[out] prob The probability of sampling the returned index.
	 * @return The index of the sampled discrete bin.
	 */
	uint32_t sample(Float u, Float &prob) const;

	/**
	 * @brief Draws a sample and linearly interpolates within the chosen interval.
	 * @param u A uniformly distributed random number in [0, 1).
	 * @param[out] pdf The probability density function (PDF) value at the sampled point.
	 * @param[out] index The index of the discrete bin chosen.
	 * @return A continuous sample value in [0, 1].
	 */
	Float sampleLinear(Float u, Float &pdf, uint32_t &index) const;

	/**
	 * @brief Returns the PDF for a given discrete bin index.
	 * @param index The index of the bin.
	 * @return The probability density function (PDF) value.
	 */
	Float pdf(uint32_t index) const;

	/**
	 * @brief Returns the PDF for a continuous position 'x' in [0, 1].
	 * @param x The continuous position.
	 * @return The probability density function (PDF) value.
	 */
	Float pdf(Float x) const;

	/// Cumulative distribution function (CDF).
	std::vector<Float> m_cdf;
	/// The integral of the original function (unnormalized).
	Float m_acc;
};

/**
 * @brief Implements a 2D discrete probability distribution for importance sampling.
 *
 * This class takes a 2D array of function values (e.g., an image) and allows
 * for efficient sampling, where samples are drawn proportionally to the values.
 *
 * It uses the method of conditional probabilities, decomposing the 2D sampling
 * problem into two 1D sampling problems: one for the marginal distribution and
 * one for the conditional distributions.
 */
struct DistribTable2D {
	DistribTable2D() = default;

	/**
	 * @brief Constructs a 2D distribution from a 2D array of weights.
	 * @param f A pointer to a 2D array of 'nx' * 'ny' values, stored row-major.
	 * @param nx The width of the 2D array.
	 * @param ny The height of the 2D array.
	 */
	DistribTable2D(const Float *f, uint32_t nx, uint32_t ny);

	/**
	 * @brief Draws a continuous 2D sample from the distribution.
	 * @param u A uniformly distributed 2D random vector with components in [0, 1).
	 * @param[out] pdf The probability density function (PDF) value at the sampled point.
	 * @return A continuous 2D sample vector with components in [0, 1].
	 */
	Vector2 sampleLinear(const Vector2 &u, Float &pdf) const;

	/**
	 * @brief Returns the PDF for a given discrete coordinate (x, y).
	 * @param x The discrete x-coordinate.
	 * @param y The discrete y-coordinate.
	 * @return The probability density function (PDF) value.
	 */
	Float pdf(uint32_t x, uint32_t y) const;

	/**
	 * @brief Returns the PDF for a continuous position 'p' in [0, 1]^2.
	 * @param p The continuous 2D position.
	 * @return The probability density function (PDF) value.
	 */
	Float pdf(const Vector2 &p) const;

	/// Conditional distributions p(x|y) for each row y.
	std::vector<DistribTable> m_cond;
	/// Marginal distribution p(y).
	DistribTable m_margin;
};

MTS_NAMESPACE_END

#endif // __DISTRIB_H