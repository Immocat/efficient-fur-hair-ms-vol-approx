/*
	This file is part of Mitsuba, a physically based rendering system.

	Copyright (c) 2007-2014 by Wenzel Jakob and others.

	Mitsuba is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License Version 3
	as published by the Free Software Foundation.

	Mitsuba is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#if !defined(__MITSUBA_CORE_SHVECTOR_H_)
#define __MITSUBA_CORE_SHVECTOR_H_

#include <Eigen/Core>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/mitsuba.h>
#include <mitsuba/render/sampler.h>
MTS_NAMESPACE_BEGIN

using DynamicMatrix = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
using DynamicVector = Eigen::Matrix<Float, Eigen::Dynamic, 1>;

/* Precompute normalization coefficients for the first 10 bands */
#define SH_NORMTBL_SIZE 10

struct SHVector;

/**
 * \brief Stores the diagonal blocks of a spherical harmonic
 * rotation matrix
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE SHRotation {
	typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> Matrix;

	std::vector<Matrix> blocks;

	/// Construct a new rotation storage for the given number of bands
	inline SHRotation(int bands)
		: blocks(bands) {
		for (int i = 0; i < bands; ++i) {
			int dim = 2 * i + 1;
			blocks[i] = Matrix(dim, dim);
		}
	}

	/**
	 * \brief Transform a coefficient vector and store the result into
	 * the given target vector.
	 *
	 * The source and target must have the same number of bands.
	 */
	void operator()(const SHVector &source, SHVector &target) const;
};

/**
 * \brief Stores a truncated real spherical harmonics representation of
 * an L2-integrable function.
 *
 * Also provides some other useful functionality, such as evaluation,
 * projection and rotation.
 *
 * The Mathematica equivalent of the basis functions implemented here is:
 *
 * \code
 * SphericalHarmonicQ[l_, m_, \[Theta]_, \[Phi]_] :=
 *   Piecewise[{
 *      {SphericalHarmonicY[l, m, \[Theta], \[Phi]], m == 0},
 *      {Sqrt[2]*Re[SphericalHarmonicY[l, m, \[Theta], \[Phi]]], m > 0},
 *      {Sqrt[2]*Im[SphericalHarmonicY[l, -m, \[Theta], \[Phi]]], m < 0}
 *  }]
 * \endcode
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct MTS_EXPORT_CORE SHVector {
  public:
	/// Construct an invalid SH vector
	inline SHVector()
		: m_bands(0) {
	}

	/// Construct a new SH vector (initialized to zero)
	inline SHVector(int bands)
		: m_bands(bands), m_coeffs(bands * bands) {
		clear();
	}

	/// Unserialize a SH vector to a binary data stream
	SHVector(Stream *stream);

	/// Copy constructor
	inline SHVector(const SHVector &v)
		: m_bands(v.m_bands),
		  m_coeffs(v.m_coeffs) {
	}

	/// Return the number of stored SH coefficient bands
	inline int getBands() const {
		return m_bands;
	}

	/// Serialize a SH vector to a binary data stream
	void serialize(Stream *stream) const;

	/// Get the energy per band
	inline Float energy(int band) const {
		Float result = 0;
		for (int m = -band; m <= band; ++m)
			result += std::abs(operator()(band, m));
		return result;
	}

	/// Assignment
	inline SHVector &operator=(const SHVector &v) {
		m_bands = v.m_bands;
		m_coeffs = v.m_coeffs;
		return *this;
	}

	/// Set all coefficients to zero
	inline void clear() {
		m_coeffs.setZero();
	}

	/// Component-wise addition
	inline SHVector &operator+=(const SHVector &v) {
		if (v.m_coeffs.size() > m_coeffs.size()) {
			ptrdiff_t oldSize = m_coeffs.size();
			m_coeffs.conservativeResize(v.m_coeffs.size());
			m_coeffs.tail(m_coeffs.size() - oldSize).setZero();
			m_bands = v.m_bands;
		}

		ptrdiff_t commonSize = std::min(m_coeffs.size(), v.m_coeffs.size());

		if (commonSize > 0) {
			m_coeffs.head(commonSize) += v.m_coeffs.head(commonSize);
		}

		return *this;
	}

	/// Component-wise addition
	inline SHVector operator+(const SHVector &v) const {
		SHVector vec(std::max(m_bands, v.m_bands));
		if (m_bands > v.m_bands) {
			vec.m_coeffs = m_coeffs;
			vec.m_coeffs.head(v.m_coeffs.size()) += v.m_coeffs;
		} else {
			vec.m_coeffs = v.m_coeffs;
			vec.m_coeffs.head(m_coeffs.size()) += m_coeffs;
		}
		return vec;
	}

	/// Component-wise subtraction
	inline SHVector &operator-=(const SHVector &v) {
		ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
		if (extendBy > 0) {
			m_coeffs.conservativeResize(v.m_coeffs.rows());
			m_coeffs.tail(extendBy).setZero();
			m_bands = v.m_bands;
		}
		m_coeffs.head(m_coeffs.size()) -= v.m_coeffs.head(m_coeffs.size());
		return *this;
	}

	/// Component-wise subtraction
	inline SHVector operator-(const SHVector &v) const {
		SHVector vec(std::max(m_bands, v.m_bands));
		if (m_bands > v.m_bands) {
			vec.m_coeffs = m_coeffs;
			vec.m_coeffs.head(v.m_coeffs.size()) -= v.m_coeffs;
		} else {
			vec.m_coeffs = -v.m_coeffs;
			vec.m_coeffs.head(m_coeffs.size()) += m_coeffs;
		}
		return vec;
	}

	/// Add a scalar multiple of another vector
	inline SHVector &madd(Float f, const SHVector &v) {
		ptrdiff_t extendBy = v.m_coeffs.size() - m_coeffs.size();
		if (extendBy > 0) {
			m_coeffs.conservativeResize(v.m_coeffs.rows());
			m_coeffs.tail(extendBy).setZero();
			m_bands = v.m_bands;
		}
		m_coeffs.head(m_coeffs.size()) += v.m_coeffs.head(m_coeffs.size()) * f;

		return *this;
	}

	/// Scalar multiplication
	inline SHVector &operator*=(Float f) {
		m_coeffs *= f;
		return *this;
	}

	/// Scalar multiplication
	inline SHVector operator*(Float f) const {
		SHVector vec(m_bands);
		vec.m_coeffs = m_coeffs * f;
		return vec;
	}

	/// Scalar division
	inline SHVector &operator/=(Float f) {
		m_coeffs *= (Float)1 / f;
		return *this;
	}

	/// Scalar division
	inline SHVector operator/(Float f) const {
		SHVector vec(m_bands);
		vec.m_coeffs = m_coeffs * (1 / f);
		return vec;
	}

	/// Negation operator
	inline SHVector operator-() const {
		SHVector vec(m_bands);
		vec.m_coeffs = -m_coeffs;
		return vec;
	}

	/// Access coefficient m (in {-l, ..., l}) on band l
	inline Float &operator()(int l, int m) {
		return m_coeffs[l * (l + 1) + m];
	}

	/// Access coefficient m (in {-l, ..., l}) on band l
	inline const Float &operator()(int l, int m) const {
		return m_coeffs[l * (l + 1) + m];
	}

	/// Evaluate for a direction given in spherical coordinates
	Float eval(Float theta, Float phi) const;

	/// Evaluate for a direction given in Cartesian coordinates
	Float eval(const Vector &v) const;

	/**
	 * \brief Evaluate for a direction given in spherical coordinates.
	 *
	 * This function is much faster but only works for azimuthally
	 * invariant functions
	 */
	Float evalAzimuthallyInvariant(Float theta, Float phi) const;

	/**
	 * \brief Evaluate for a direction given in cartesian coordinates.
	 *
	 * This function is much faster but only works for azimuthally
	 * invariant functions
	 */
	Float evalAzimuthallyInvariant(const Vector &v) const;

	/// Check if this function is azumuthally invariant
	bool isAzimuthallyInvariant() const;

	/// Equality comparison operator
	inline bool operator==(const SHVector &v) const {
		return m_bands == v.m_bands && m_coeffs == v.m_coeffs;
	}

	/// Equality comparison operator
	inline bool operator!=(const SHVector &v) const {
		return !operator==(v);
	}

	/// Dot product
	inline friend Float dot(const SHVector &v1, const SHVector &v2);

	/// Normalize so that the represented function becomes a valid distribution
	void normalize();

	/// Compute the second spherical moment (analytic)
	Matrix3x3 mu2() const;

	/// Brute-force search for the minimum value over the sphere
	Float findMinimum(int res) const;

	/// Add a constant value
	void addOffset(Float value);

	/**
	 * \brief Convolve the SH representation with the supplied kernel.
	 *
	 * Based on the Funk-Hecke theorem -- the kernel must be rotationally
	 * symmetric around the Z-axis.
	 */
	void convolve(const SHVector &kernel);

	/// Project the given function onto a SH basis (using a 2D composite Simpson's rule)
	template <typename Functor>
	void project(const Functor &f, int res = 32) {
		SAssert(res % 2 == 0);
		/* Nested composite Simpson's rule */
		Float hExt = M_PI / res,
			  hInt = (2 * M_PI) / (res * 2);

		for (int l = 0; l < m_bands; ++l)
			for (int m = -l; m <= l; ++m)
				operator()(l, m) = 0;

		Float *sinPhi = (Float *)alloca(sizeof(Float) * m_bands),
			  *cosPhi = (Float *)alloca(sizeof(Float) * m_bands);

		for (int i = 0; i <= res; ++i) {
			Float theta = hExt * i, cosTheta = std::cos(theta);
			Float weightExt = (i & 1) ? 4.0f : 2.0f;
			if (i == 0 || i == res)
				weightExt = 1.0f;

			for (int j = 0; j <= res * 2; ++j) {
				Float phi = hInt * j;
				Float weightInt = (j & 1) ? 4.0f : 2.0f;
				if (j == 0 || j == 2 * res)
					weightInt = 1.0f;

				for (int m = 0; m < m_bands; ++m) {
					sinPhi[m] = std::sin((m + 1) * phi);
					cosPhi[m] = std::cos((m + 1) * phi);
				}

				Float value = f(sphericalDirection(theta, phi)) * std::sin(theta) * weightExt * weightInt;

				for (int l = 0; l < m_bands; ++l) {
					for (int m = 1; m <= l; ++m) {
						Float L = legendreP(l, m, cosTheta) * normalization(l, m);
						operator()(l, -m) += value * SQRT_TWO * sinPhi[m - 1] * L;
						operator()(l, m) += value * SQRT_TWO * cosPhi[m - 1] * L;
					}

					operator()(l, 0) += value * legendreP(l, 0, cosTheta) * normalization(l, 0);
				}
			}
		}

		for (int l = 0; l < m_bands; ++l)
			for (int m = -l; m <= l; ++m)
				operator()(l, m) *= hExt * hInt / 9;
	}

	/// Compute the relative L2 error
	template <typename Functor>
	Float l2Error(const Functor &f, int res = 32) const {
		SAssert(res % 2 == 0);
		/* Nested composite Simpson's rule */
		Float hExt = M_PI / res,
			  hInt = (2 * M_PI) / (res * 2);
		Float error = 0.0f, denom = 0.0f;

		for (int i = 0; i <= res; ++i) {
			Float theta = hExt * i;
			Float weightExt = (i & 1) ? 4.0f : 2.0f;
			if (i == 0 || i == res)
				weightExt = 1.0f;

			for (int j = 0; j <= res * 2; ++j) {
				Float phi = hInt * j;
				Float weightInt = (j & 1) ? 4.0f : 2.0f;
				if (j == 0 || j == 2 * res)
					weightInt = 1.0f;

				Float value1 = f(sphericalDirection(theta, phi));
				Float value2 = eval(theta, phi);
				Float diff = value1 - value2;
				Float weight = std::sin(theta) * weightInt * weightExt;

				error += diff * diff * weight;
				denom += value1 * value1 * weight;
			}
		}

		return error / denom;
	}

	/**
	 * \brief Transform the coefficients by a general matrix M.
	 * Result = M * this.
	 * M must be (bands*bands) x (bands*bands).
	 */
	inline void applyMatrix(const Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> &M, SHVector &target) const {
		SAssert(M.cols() == m_coeffs.rows());

		// Resize target if needed
		if (target.m_bands != m_bands) {
			target.m_bands = m_bands;
			target.m_coeffs.resize(m_coeffs.rows());
		}

		// Direct Eigen multiplication
		target.m_coeffs = M * m_coeffs;
	}

	/// Turn into a string representation
	std::string toString() const;

	/// Return a normalization coefficient
	inline static Float normalization(int l, int m) {
		if (l < SH_NORMTBL_SIZE)
			return m_normalization[l * (l + 1) / 2 + m];
		else
			return computeNormalization(l, m);
	}

	/**
	 * \brief Recursively computes rotation matrices for each band of SH coefficients.
	 *
	 * Based on 'Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion'
	 * by Ivanic and Ruedenberg. The implemented tables follow the notation in
	 * 'Spherical Harmonic Lighting: The Gritty Details' by Robin Green.
	 */
	static void rotation(const Transform &t, SHRotation &rot);

	/// Precomputes normalization coefficients for the first few bands
	static void staticInitialization();

	/// Free the memory taken up by staticInitialization()
	static void staticShutdown();

  protected:
	/// Helper function for rotation() -- computes a diagonal block based on the previous level
	static void rotationBlock(const SHRotation::Matrix &M1, const SHRotation::Matrix &Mp, SHRotation::Matrix &Mn);

	/// Compute a normalization coefficient
	static Float computeNormalization(int l, int m);

  private:
	int m_bands;
	Eigen::Matrix<Float, Eigen::Dynamic, 1> m_coeffs;
	static Float *m_normalization;
};

inline Float dot(const SHVector &v1, const SHVector &v2) {
	const size_t size = std::min(v1.m_coeffs.size(), v2.m_coeffs.size());
	return v1.m_coeffs.head(size).dot(v2.m_coeffs.head(size));
}

/**
 * \brief Implementation of 'Importance Sampling Spherical Harmonics'
 * by W. Jarsz, N. Carr and H. W. Jensen (EUROGRAPHICS 2009)
 *
 * \ingroup libcore
 * \ingroup libpython
 */
class MTS_EXPORT_CORE SHSampler : public Object {
  public:
	/**
	 * \brief Precompute a spherical harmonics sampler object
	 *
	 * \param bands Number of SH coefficient bands to support
	 * \param depth Number of recursive sample warping steps.
	 */
	SHSampler(int bands, int depth);

	/**
	 * \brief Warp a uniform sample in [0,1]^2 to one that is
	 * approximately proportional to the specified function.
	 *
	 * The resulting sample will have spherical coordinates
	 * [0,pi]x[0,2pi] and its actual PDF (which might be
	 * slightly different from the function evaluated at the
	 * sample, even if $f$ is a distribution) will be returned.
	 */
	Float warp(const SHVector &f, Point2 &sample) const;

	/// Return information on the size of the precomputed tables
	std::string toString() const;

	MTS_DECLARE_CLASS()
  protected:
	/// Virtual destructor
	virtual ~SHSampler();

	/* Index into the assoc. legendre polynomial table */
	inline int I(int l, int m) const {
		return l * (l + 1) / 2 + m;
	}

	/* Index into the phi table */
	inline int P(int m) const {
		return m + m_bands;
	}

	inline Float lookupIntegral(int depth, int zBlock, int phiBlock, int l, int m) const {
		return -m_phiMap[depth][phiBlock][P(m)] * m_legendreMap[depth][zBlock][I(l, std::abs(m))];
	}

	/// Recursively compute assoc. legendre & phi integrals
	Float *legendreIntegrals(Float a, Float b);
	Float *phiIntegrals(Float a, Float b);

	/// Integrate a SH expansion over the specified mip-map region
	Float integrate(int depth, int zBlock, int phiBlock, const SHVector &f) const;

  protected:
	int m_bands;
	int m_depth;
	Float ***m_phiMap;
	Float ***m_legendreMap;
	int m_dataSize;
	Float *m_normalization;
};

struct MatrixRGB {
	DynamicMatrix r;
	DynamicMatrix g;
	DynamicMatrix b;
};

struct SHVectorRGB {
	SHVector r, g, b;

	SHVectorRGB() = default;

	inline SHVectorRGB(int bands)
		: r(bands), g(bands), b(bands) {
	}

	/// Apply a different matrix to each color channel.
	inline void applyMatrices(
		const DynamicMatrix &redMatrix,
		const DynamicMatrix &greenMatrix,
		const DynamicMatrix &blueMatrix,
		SHVectorRGB &target) const {
		r.applyMatrix(redMatrix, target.r);
		g.applyMatrix(greenMatrix, target.g);
		b.applyMatrix(blueMatrix, target.b);
	}

	/// Apply the same matrix to all color channels.
	inline void applyMatrix(const DynamicMatrix &matrix, SHVectorRGB &target) const {
		r.applyMatrix(matrix, target.r);
		g.applyMatrix(matrix, target.g);
		b.applyMatrix(matrix, target.b);
	}

	inline SHVectorRGB &operator+=(const SHVectorRGB &value) {
		r += value.r;
		g += value.g;
		b += value.b;
		return *this;
	}

	/// Multiply each channel by a scalar weight.
	inline SHVectorRGB operator*(Float factor) const {
		SHVectorRGB result(r.getBands());
		result.r = r * factor;
		result.g = g * factor;
		result.b = b * factor;
		return result;
	}

	/// Add two RGB SH expansions component-wise.
	inline SHVectorRGB operator+(const SHVectorRGB &other) const {
		SHVectorRGB result(std::max(r.getBands(), other.r.getBands()));
		result.r = r + other.r;
		result.g = g + other.g;
		result.b = b + other.b;
		return result;
	}

	/**
	 * \brief Evaluate the RGB SH expansion for a direction.
	 *
	 * Each channel is evaluated independently. Negative values are clamped
	 * to zero to suppress ringing artifacts from the truncated SH basis.
	 */
	inline Spectrum eval(const Vector &dir) const {
		const Float redValue = r.eval(dir);
		const Float greenValue = g.eval(dir);
		const Float blueValue = b.eval(dir);

		Spectrum result;
		result.fromLinearRGB(
			std::max((Float) 0.0f, redValue),
			std::max((Float) 0.0f, greenValue),
			std::max((Float) 0.0f, blueValue));
		return result;
	}

	/**
	 * \brief Project a directional RGB sample into the SH basis.
	 *
	 * The temporary vector is used as a reusable scratch buffer containing
	 * the evaluated basis values for the queried direction.
	 */
	inline void project(const Vector &dir, const Spectrum &energy, SHVector &basisValues) {
		const int bands = basisValues.getBands();
		if (getBands() != bands)
			*this = SHVectorRGB(bands);
		const Float cosTheta = dir.z;
		Float phi = std::atan2(dir.y, dir.x);
		if (phi < 0)
			phi += 2 * (Float) M_PI;

		for (int l = 0; l < bands; ++l) {
			for (int m = -l; m <= l; ++m) {
				const int absM = std::abs(m);
				const Float normalizationFactor = SHVector::normalization(l, absM);
				const Float legendreValue = legendreP(l, absM, cosTheta);
				Float basisValue;
				if (m == 0)
					basisValue = normalizationFactor * legendreValue;
				else if (m > 0)
					basisValue = SQRT_TWO * normalizationFactor * legendreValue * std::cos(m * phi);
				else
					basisValue = SQRT_TWO * normalizationFactor * legendreValue * std::sin(absM * phi);
				basisValues(l, m) = basisValue;
				r(l, m) = energy[0] * basisValue;
				g(l, m) = energy[1] * basisValue;
				b(l, m) = energy[2] * basisValue;
			}
		}
	}

	inline int getBands() const {
		return r.getBands();
	}
};

/**
 * \brief Rotates an RGB SH vector using a precomputed SH rotation matrix.
 *
 * \param rot    The precomputed SHRotation object.
 * \param source The source SHVectorRGB (e.g., World Space).
 * \param target The target SHVectorRGB (e.g., Local Space).
 */
void MTS_EXPORT_CORE RotateSHRGB(const SHRotation &rot, const SHVectorRGB &source, SHVectorRGB &target);

/**
 * \brief Creates a Spherical Harmonics representation of a Clamped Cosine Lobe
 * oriented in a specific direction.
 *
 * This implementation uses the EXACT analytical Zonal Harmonics coefficients
 * for the function f(theta) = max(cos(theta), 0).
 *
 *
 * \param direction  Target direction. It will be normalized internally.
 * \param peakRadiance  The scalar radiance to scale the lobe by.
 * \param bands      Number of SH bands (e.g., 3 bands = 9 coefficients).
 * \param outSH      The output SHVector to store the result.
 */
void MTS_EXPORT_CORE createScaledCosineLobe(const Vector &direction, Float peakRadiance, int bands, SHVector &outSH);

// Create Cosine Lobe according to LPV paper
SHVector MTS_EXPORT_CORE createLPVCosineLobe(const Vector &direction);

SHVector MTS_EXPORT_CORE createDirectionRay(const Vector3 &direction, int bands);

SHVector MTS_EXPORT_CORE createDirectionLobe(const Vector3 &direction, int bands);
/**
 * \brief Evaluates a single Spherical Harmonic basis function Y_lm(dir).
 *
 * \param index The linear SH index (index = l * (l + 1) + m)
 * \param dir   The normalized direction vector in Cartesian coordinates.
 * \return      The value of the SH basis function.
 */
Float MTS_EXPORT_CORE EvalSH(int index, const Vector &dir);
/**
 * \brief Efficiently evaluates all SH basis functions up to the specified number of bands.
 *
 * This function is optimized for Monte Carlo loops. It avoids re-calculating
 * trigonometric functions for the same direction.
 *
 * \param bands     The number of SH bands (e.g., 3 bands = 9 coeffs).
 * \param dir       The normalized direction vector.
 * \param results   Output vector. Will be resized automatically.
 *                  Passing by reference allows memory reuse across MC iterations.
 */
void MTS_EXPORT_CORE EvalSH_All(int bands, const Vector &dir, std::vector<Float> &results);

void MTS_EXPORT_CORE EvalSH_All_Eigen(int bands, const Vector &dir, DynamicVector &results);

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_SHVECTOR_H_ */
