#include "mitsuba/core/distrib.h"
#include <numeric> // For std::accumulate if needed

MTS_NAMESPACE_BEGIN

DistribTable::DistribTable(const Float *f, uint32_t n)
	: m_cdf(n + 1) {
	m_cdf[0] = 0;
	Float invn = 1.0f / (Float)n;
	for (uint32_t i = 1; i < n + 1; ++i) {
		// Build the unnormalized CDF by summing up function values.
		m_cdf[i] = m_cdf[i - 1] + f[i - 1] * invn;
	}

	// The accumulator 'm_acc' stores the integral of the function.
	m_acc = m_cdf[n];

	if (m_acc == 0.0f) {
		// If the integral is zero, fall back to a uniform distribution.
		for (uint32_t i = 0; i < n + 1; ++i) {
			m_cdf[i] = (Float)i / (Float)n;
		}
	} else {
		// Normalize the CDF so that its last element is 1.
		Float invAcc = 1.0f / m_acc;
		for (uint32_t i = 0; i < n + 1; ++i) {
			m_cdf[i] *= invAcc;
		}
	}
	// Ensure the last element is exactly 1.0 to avoid numerical precision issues.
	m_cdf[n] = 1.0f;
}

uint32_t DistribTable::sample(Float u, Float &prob) const {
	u = std::clamp(u, 0.0f, std::nextafter(1.0f, 0.0f));
	assert(u >= 0.0f && u < 1.0f);

	// Find the first CDF entry that is greater than 'u' using binary search.
	auto it = std::upper_bound(m_cdf.begin(), m_cdf.end(), u);
	uint32_t index = (uint32_t)std::distance(m_cdf.begin(), std::prev(it));

	// The probability is the size of the interval in the CDF.
	prob = m_cdf[index + 1] - m_cdf[index];
	return index;
}

Float DistribTable::sampleLinear(Float u, Float &pdf, uint32_t &index) const {
	u = std::clamp(u, 0.0f, std::nextafter(1.0f, 0.0f));
	assert(u >= 0.0f && u < 1.0f);

	// Find the discrete bin for the random number 'u'.
	auto it = std::upper_bound(m_cdf.begin(), m_cdf.end(), u);
	index = (uint32_t)std::distance(m_cdf.begin(), std::prev(it));

	Float n = (Float)(m_cdf.size() - 1);
	// PDF = Probability / IntervalWidth = (cdf[i+1] - cdf[i]) / (1/n)
	pdf = (m_cdf[index + 1] - m_cdf[index]) * n;

	// Linearly interpolate within the selected interval.
	Float du = (u - m_cdf[index]) / (m_cdf[index + 1] - m_cdf[index]);
	return (index + du) / n;
}

Float DistribTable::pdf(uint32_t index) const {
	uint32_t n = (uint32_t)(m_cdf.size() - 1);
	assert(index < n);
	// PDF = Probability / IntervalWidth
	return (m_cdf[index + 1] - m_cdf[index]) * (Float)n;
}

Float DistribTable::pdf(Float x) const {
	assert(x >= 0.0 && x <= 1.0f);
	int n = (int)(m_cdf.size() - 1);
	uint32_t idx = (uint32_t)std::clamp((int)std::floor(x * n), 0, n - 1);
	return pdf(idx);
}

// =======================================================================

DistribTable2D::DistribTable2D(const Float *f, uint32_t nx, uint32_t ny) {
	m_cond.resize(ny);
	// 1. Create a 1D conditional distribution p(x|y) for each row y.
	for (uint32_t y = 0; y < ny; ++y) {
		m_cond[y] = DistribTable(f + (y * nx), nx);
	}

	// 2. Create the 1D marginal distribution p(y). The weight for each y
	//    is the integral of the corresponding row, stored in m_cond[y].m_acc.
	std::vector<Float> marginPmf(ny);
	for (uint32_t y = 0; y < ny; ++y) {
		marginPmf[y] = m_cond[y].m_acc;
	}
	m_margin = DistribTable(marginPmf.data(), ny);
}

Vector2 DistribTable2D::sampleLinear(const Vector2 &u, Float &pdf) const {
	Float pdfMargin, pdfCond;
	uint32_t uIndex, vIndex;

	// Sample the marginal distribution p(y) to select a row.
	Float vSample = m_margin.sampleLinear(u[1], pdfMargin, vIndex);
	// Given the row, sample the conditional distribution p(x|y).
	Float uSample = m_cond[vIndex].sampleLinear(u[0], pdfCond, uIndex);

	// The joint PDF is the product: pdf(x,y) = pdf(x|y) * pdf(y)
	pdf = pdfMargin * pdfCond;
	return Vector2(uSample, vSample);
}

Float DistribTable2D::pdf(uint32_t x, uint32_t y) const {
	// pdf(x,y) = pdf(x|y) * pdf(y)
	return m_margin.pdf(y) * m_cond[y].pdf(x);
}

Float DistribTable2D::pdf(const Vector2 &p) const {
	int x = std::clamp(int(p[0] * (m_cond[0].m_cdf.size() - 1)), 0, (int)m_cond[0].m_cdf.size() - 2);
	int y = std::clamp(int(p[1] * m_cond.size()), 0, (int)m_cond.size() - 1);
	return pdf(x, y);
}

MTS_NAMESPACE_END