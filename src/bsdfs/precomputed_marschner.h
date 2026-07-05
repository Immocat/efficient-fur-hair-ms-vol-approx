#pragma once
#if !defined(__PRECOMPUTED_MARSCHNER_H)
#define __PRECOMPUTED_MARSCHNER_H
#include <mitsuba/core/distrib.h>
#include <mitsuba/mitsuba.h>
MTS_NAMESPACE_BEGIN
// #pragma optimize(off, "")
class InterpolatedDistribution1D : public Object {
	int _size;
	int _numDistributions;

	std::vector<Float> _pdfs;
	std::vector<Float> _cdfs;
	std::vector<Float> _sums;

	Float cdfs(int x, int distribution) const {
		return _cdfs[x + distribution * (_size + 1)];
	}

	Float pdfs(int x, int distribution) const {
		return _pdfs[x + distribution * _size];
	}

	Float &cdfs(int x, int distribution) {
		return _cdfs[x + distribution * (_size + 1)];
	}

	Float &pdfs(int x, int distribution) {
		return _pdfs[x + distribution * _size];
	}

  public:
	InterpolatedDistribution1D(std::vector<Float> weights, int size, int numDistributions)
		: _size(size), _numDistributions(numDistributions), _pdfs(std::move(weights)),
		  _cdfs((size + 1) * numDistributions), _sums(numDistributions) {
		for (int dist = 0; dist < _numDistributions; ++dist) {
			cdfs(0, dist) = 0.0f;
			for (int x = 0; x < _size; ++x)
				cdfs(x + 1, dist) = pdfs(x, dist) + cdfs(x, dist);

			_sums[dist] = cdfs(_size, dist);

			if (_sums[dist] < 1e-4f) {
				// Revert to uniform sampling for near-degenerate distributions
				Float ratio = 1.0f / _size;
				for (int x = 0; x < _size; ++x) {
					pdfs(x, dist) = ratio;
					cdfs(x, dist) = x * ratio;
				}
			} else {
				Float scale = 1.0f / _sums[dist];
				for (int x = 0; x < _size; ++x) {
					pdfs(x, dist) *= scale;
					cdfs(x, dist) *= scale;
				}
			}
			cdfs(_size, dist) = 1.0f;
		}
	}

	void warp(Float distribution, Float &u, int &x) const {
		int d0 = math::clamp(int(distribution), 0, _numDistributions - 1);
		int d1 = std::min(d0 + 1, _numDistributions - 1);
		Float v = math::clamp(distribution - d0, 0.0f, 1.0f);

		int lower = 0, upper = _size;
		Float lowerU = 0.0f, upperU = 1.0f;

		while (upper - lower != 1) {
			int midpoint = (upper + lower) / 2;
			Float midpointU = cdfs(midpoint, d0) * (1.0f - v) + cdfs(midpoint, d1) * v;
			if (midpointU < u) {
				lower = midpoint;
				lowerU = midpointU;
			} else {
				upper = midpoint;
				upperU = midpointU;
			}
		}

		x = lower;
		u = math::clamp((u - lowerU) / (upperU - lowerU), 0.0f, 1.0f);
	}

	Float pdf(Float distribution, int x) const {
		int d0 = math::clamp(int(distribution), 0, _numDistributions - 1);
		int d1 = std::min(d0 + 1, _numDistributions - 1);
		Float v = math::clamp(distribution - d0, 0.0f, 1.0f);

		return pdfs(x, d0) * (1.0f - v) + pdfs(x, d1) * v;
	}

	Float sum(Float distribution) const {
		int d0 = math::clamp(int(distribution), 0, _numDistributions - 1);
		int d1 = std::min(d0 + 1, _numDistributions - 1);
		Float v = math::clamp(distribution - d0, 0.0f, 1.0f);

		return _sums[d0] * (1.0f - v) + _sums[d1] * v;
	}
	MTS_DECLARE_CLASS()
};
MTS_IMPLEMENT_CLASS(InterpolatedDistribution1D, false, Object)

class PrecomputedAzimuthalLobe : public Object {
  public:
	static const int AzimuthalResolution = 64;

  public:
	std::vector<Vector3> _table;
	std::unique_ptr<InterpolatedDistribution1D> _sampler;

  public:
	PrecomputedAzimuthalLobe(std::vector<Vector3> table)
		: _table(std::move(table)) {
		const int Size = AzimuthalResolution;

		std::vector<float> weights(Size * Size);
		for (int i = 0; i < Size * Size; ++i)
			weights[i] = _table[i].length();

		// Dilate weights slightly to stay conservative
		for (int y = 0; y < Size; ++y) {
			for (int x = 0; x < Size - 1; ++x)
				weights[x + y * Size] = std::max(weights[x + y * Size], weights[x + 1 + y * Size]);
			for (int x = Size - 1; x > 0; --x)
				weights[x + y * Size] = std::max(weights[x + y * Size], weights[x - 1 + y * Size]);
		}
		for (int x = 0; x < Size; ++x) {
			for (int y = 0; y < Size - 1; ++y)
				weights[x + y * Size] = std::max(weights[x + y * Size], weights[x + (y + 1) * Size]);
			for (int y = Size - 1; y > 0; --y)
				weights[x + y * Size] = std::max(weights[x + y * Size], weights[x + (y - 1) * Size]);
		}

		_sampler.reset(new InterpolatedDistribution1D(std::move(weights), Size, Size));
	}

	void sample(float cosThetaD, float xi, float &phi, float &pdf) const {
		float v = (AzimuthalResolution - 1) * cosThetaD;

		int x;
		_sampler->warp(v, xi, x);
		if (isinf(xi) || isnan(xi)) {
			SLog(EWarn, "PrecomputedAzimuthalLobe::sample(): non-finite warped sample; clamping to zero.");
			xi = 0.f;
		}
		phi = 2 * M_PI * (x + xi) * (1.0f / AzimuthalResolution);
		if (isinf(phi) || isnan(phi)) {
			SLog(EWarn, "PrecomputedAzimuthalLobe::sample(): non-finite azimuth; returning zero PDF.");
			phi = 0.f;
			pdf = 0.f;
			return;
		}

		pdf = _sampler->pdf(v, x) * float(AzimuthalResolution * 0.5f * INV_PI);
	}

	Vector3 eval(float phi, float cosThetaD) const {
		float u = (AzimuthalResolution - 1) * phi * 0.5f * INV_PI;
		float v = (AzimuthalResolution - 1) * cosThetaD;
		int x0 = math::clamp(int(u), 0, AzimuthalResolution - 2);
		int y0 = math::clamp(int(v), 0, AzimuthalResolution - 2);
		int x1 = x0 + 1;
		int y1 = y0 + 1;
		u = math::clamp(u - x0, 0.0f, 1.0f);
		v = math::clamp(v - y0, 0.0f, 1.0f);

		return (_table[x0 + y0 * AzimuthalResolution] * (1.0f - u) + _table[x1 + y0 * AzimuthalResolution] * u) *
				   (1.0f - v) +
			   (_table[x0 + y1 * AzimuthalResolution] * (1.0f - u) + _table[x1 + y1 * AzimuthalResolution] * u) * v;
	}

	float pdf(float phi, float cosThetaD) const {
		float u = (AzimuthalResolution - 1) * phi * 0.5f * INV_PI;
		float v = (AzimuthalResolution - 1) * cosThetaD;

		int x = int(u);
		if (x < 0 || x >= AzimuthalResolution) {
			Log(EError, "Phi overflow. Phi: %f, x: %d", phi, x);
		}
		x = math::clamp(x, 0, AzimuthalResolution - 1);

		return _sampler->pdf(v, x) * float(AzimuthalResolution * 0.5f * INV_PI);
	}

	float weight(float cosThetaD) const {
		float v = (AzimuthalResolution - 1) * cosThetaD;
		return _sampler->sum(v) * (2 * M_PI / AzimuthalResolution);
	}
	MTS_DECLARE_CLASS()
};
MTS_IMPLEMENT_CLASS(PrecomputedAzimuthalLobe, false, Object)

// for reference of practical bysdf
class PrecomputedLongitudinalLobe {
  public:
	static const int AzimuthalResolution = 256;

  private:
	std::vector<Vector3> _table;
	DistribTable _sampler;

  public:
	PrecomputedLongitudinalLobe(std::vector<Vector3> table)
		: _table(std::move(table)) {
		const uint32_t Size = 256;

		std::vector<float> weights(Size);
		for (int i = 0; i < Size; ++i)
			weights[i] = _table[i].length();

		_sampler = DistribTable(weights.data(), Size);
	}

	void sample(float u, float &theta_i, float &pdf) const {
		uint32_t index;
		float pos = _sampler.sampleLinear(u, pdf, index);
		theta_i = pos * M_PI - 0.5 * M_PI;
		pdf /= M_PI;
	}

	Vector3 eval(float theta_i) const {
		uint32_t n = (uint32_t)(_sampler.m_cdf.size() - 1);
		uint32_t idx = uint32_t((theta_i + 0.5 * M_PI) / M_PI * n);
		return _table[idx];
	}

	float pdf(float theta_i) const {
		uint32_t n = (uint32_t)(_sampler.m_cdf.size() - 1);
		uint32_t idx = uint32_t((theta_i + 0.5 * M_PI) / M_PI * n);
		return _sampler.pdf(idx) / M_PI;
	}
	MTS_DECLARE_CLASS()
};
// #pragma optimize(on, "")
MTS_IMPLEMENT_CLASS(PrecomputedLongitudinalLobe, false, Object)

MTS_NAMESPACE_END
#endif
