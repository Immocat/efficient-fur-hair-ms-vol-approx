#pragma once
#if !defined(__MARSCHNER_H)
#define __MARSCHNER_H

#include "precomputed_marschner.h"
#include <array>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/core/shvector.h>
MTS_NAMESPACE_BEGIN

#ifdef _WIN32
#ifdef DEON_EXPORTS
#define DEON_API __declspec(dllexport)
#else
#define DEON_API __declspec(dllimport)
#endif
#else
#define DEON_API __attribute__((visibility("default")))
#endif

/**
 * @brief Implements a hair BSDF based on the d'Eon model.
 *
 * This class implements the hair scattering model described in "An Energy-Conserving
 * Hair Reflectance Model" by d'Eon et al. It incorporates the key improvements from the d'Eon paper, namely:
 *
 * 1. An energy-conserving longitudinal scattering function M_p (Equation 7).
 * 2. A unified azimuthal scattering function N_p based on numerical integration (Equation 10).
 *
 * The model decomposes hair scattering into three primary lobes:
 * - R: Direct reflection from the cuticle surface (p=0).
 * - TT: Transmission through the fiber (p=1).
 * - TRT: Transmission, internal reflection, then transmission (p=2).
 *
 * For efficiency, this implementation relies on precomputed tables for the azimuthal
 * scattering component N_p, as recommended by the paper.
 */
class DEON_API dEonHairBSDF : public BSDF {
  public:
	dEonHairBSDF(const Properties &props);
	dEonHairBSDF(Stream *stream, InstanceManager *manager);

	// Standard BSDF interface methods from Mitsuba. Implementations not shown.
	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const override;
	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const override;
	Spectrum sampleWithWorldWiAndOrientation(const Vector3 &worldWi, const Vector3 &orientation, Float &pdf, Vector3f &worldWo, Sampler *sampler) const;

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const;
	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const;
	// The same as sample, but only for debug purposes
	Vector3f test_sample(const Vector3f &wi, const Point2f &controlled_xiN, const Point2f &controlled_xiM) const;

	// wd is already in local frame
	Spectrum forwardScattering(const Scene &scene, const Float &hairThickness, const BSDFSamplingRecord &bRec, const Vector3 &wd, Float dist, const Point &safeOrigin, Vector3 &sigma_f, Vector3 &trans, bool &direct) const;
	Vector3 forwardTrans(const Scene &scene, const BSDFSamplingRecord &bRec,
						 const Vector3 &wi, const Vector3 &wd,
						 const Float hairThickness, Float dist,
						 const Point &safeOrigin,
						 Vector3 &var, bool &direct) const;
	Spectrum backwardScattering(const BSDFSamplingRecord &bRec, const Vector3 &sigmaf, const bool direct, const Vector3 &trans, Vector3 &sigmab) const;

	MatrixRGB PrecomputeHairMatrixRGB(int bands, int nSamples) const;

	void configure();

	MTS_DECLARE_CLASS()
  private:
	/// @brief Index of refraction of the hair fiber. Typically ~1.55 for human hair.
	/// Corresponds to the symbol η (eta) in the paper.
	const Float Eta = 1.55;

	/// @brief Backward scattering density factor
	const Float db = 0.7;

	/// @brief Forward scattering density factor
	const Float df = 0.7;

	// --- User-Configurable Appearance Parameters ---
	/// @brief Longitudinal shift of scattering lobes due to cuticle angle, in degrees.
	Float _scaleAngleDeg;

	/// @brief Ratio of eumelanin (brown/black) to pheomelanin (red/yellow). Range: [0, 1].
	/// Used to procedurally calculate the absorption coefficient _sigmaA.
	Float _melaninRatio;

	/// @brief Overall concentration of melanin pigments. Controls hair darkness.
	/// Also used to procedurally calculate _sigmaA.
	Float _melaninConcentration;

	/// @brief The volume absorption coefficient [m^-1], defining the hair color.
	/// This corresponds to the symbol µₐ in the paper and is crucial for calculating
	/// the attenuation term A(p,h) for internal lobes (p >= 1).
	Vector3 _sigmaA;

	/// @brief Isotropic roughness of the hair fiber. Controls the width of highlights.
	/// This value is used to derive per-lobe roughness values (β_p).
	Float _roughness;

	// --- Internal/Derived Model Parameters ---
	// These are pre-calculated in the constructor for efficiency.

	/// @brief Cuticle angle in radians, converted from _scaleAngleDeg.
	Float _scaleAngleRad;

	/// @brief Per-lobe roughness standard deviation for R, TT, and TRT lobes.
	/// Corresponds to β_p in the paper. These control the blurriness of highlights.
	Float _betaR, _betaTT, _betaTRT;

	/// @brief Per-lobe roughness variance, calculated as v_p = β_p^2.
	/// This is the direct input 'v' to the d'Eon longitudinal scattering function M_p.
	Float _vR, _vTT, _vTRT;

	/// @brief Pointers to precomputed tables for the azimuthal lobes N_p.
	/// This is a key optimization strategy from the paper to avoid costly numerical
	/// integration during rendering.
	ref<PrecomputedAzimuthalLobe> _NR, _NTT, _NTRT;

	static constexpr int NUM_BINS = 64;
	std::array<Float, NUM_BINS> alpha_b, alpha_f;
	std::array<Float, NUM_BINS> beta_b_sqr, beta_f_sqr;
	std::array<Vector3, NUM_BINS> a_b, a_f;
	std::array<Vector3, NUM_BINS> A_b, A_f, sigma_b_sqr;
	std::array<Vector3, NUM_BINS> Delta_b;

	bool generateIntegrationTable;
	bool useUniformSphereSampling = false;

	bool m_maskMode = false;
	bool m_overrideSigmaA = false;
	/**
	 * @brief Approximates the modified Bessel function of the first kind, I₀(x).
	 *
	 * Implements Equation (15) from the paper using a 10-term series expansion, which
	 * provides sufficient accuracy for rendering.
	 *
	 * @param x The input value to the function.
	 * @return The approximate value of I₀(x).
	 */
	Float I0(Float x) const;

	/**
	 * @brief Computes the natural logarithm of I₀(x).
	 *
	 * Provides a numerically stable path for large 'x' to avoid floating-point overflow
	 * that would occur when computing I₀(x) directly.
	 * @param x The input value.
	 * @return The value of log(I₀(x)).
	 */
	Float logI0(Float x) const;

	/**
	 * @brief A standard 1D Gaussian function.
	 *
	 * This is the function g(β; x) used within the 'Gaussian detector' D_p.
	 * @param beta The standard deviation of the Gaussian.
	 * @param theta The input variable.
	 * @return The value of the Gaussian PDF at theta.
	 */
	Float g(Float beta, Float theta) const;

	/**
	 * @brief Implements the "Gaussian detector" D_p from Equation (11).
	 *
	 * This function models the azimuthal blurring due to roughness. It's constructed
	 * as a sum of wrapped Gaussians to handle the periodic nature of angles.
	 * The infinite sum from the paper is approximated by a finite one.
	 *
	 * @param beta The roughness parameter (standard deviation β_p).
	 * @param phi The azimuthal angle difference.
	 * @return The value of the detector function D_p(phi).
	 */
	Float D(Float beta, Float phi) const;

	/**
	 * @brief Calculates the azimuthal exit angle for a smooth cylinder.
	 *
	 * Implements Equation (8) from the paper, which describes the geometry of light
	 * paths for different lobes.
	 *
	 * @param gammaI The arcsin of the normalized distance from the fiber center (h).
	 * @param gammaT The refracted angle corresponding to gammaI.
	 * @param p The lobe index (0 for R, 1 for TT, 2 for TRT, etc.).
	 * @return The exit angle Φ(p,h).
	 */
	Float Phi(Float gammaI, Float gammaT, int p) const;

	/**
	 * @brief The longitudinal scattering function M_p.
	 *
	 * This is a direct implementation of the paper's first key contribution,
	 * Equation (7), providing an energy-conserving model for longitudinal scattering.
	 *
	 * @param v The roughness variance, v = β².
	 * @param sinThetaI Sine of the incident longitudinal angle.
	 * @param sinThetaO Sine of the outgoing longitudinal angle.
	 * @param cosThetaI Cosine of the incident longitudinal angle.
	 * @param cosThetaO Cosine of the outgoing longitudinal angle.
	 * @return The value of M_p.
	 */
	Float M(Float v, Float sinThetaI, Float sinThetaO, Float cosThetaI, Float cosThetaO) const;

	/**
	 * @brief Importance samples the longitudinal scattering distribution M_p.
	 *
	 * This function is used by a path tracer's sampling routine to generate an
	 * outgoing direction biased towards the shape of the M_p distribution.
	 * It uses a numerically stable inversion method.
	 *
	 * @param v The roughness variance v_p.
	 * @param sinThetaI Sine of the incident longitudinal angle.
	 * @param cosThetaI Cosine of the incident longitudinal angle.
	 * @param xi1 A random number in [0, 1].
	 * @param xi2 A random number in [0, 1].
	 * @return The sine of the sampled outgoing angle, sin(θ_r).
	 */
	Float sampleM(Float v, Float sinThetaI, Float cosThetaI, Float xi1, Float xi2) const;

	/**
	 * @brief Computes the integrand for the azimuthal R-lobe (p=0), N₀.
	 * @note For reference/precomputation only. Not used at runtime.
	 *
	 * This function represents the term inside the integral of Equation (10) for p=0:
	 *   A(0,h) * D₀(φ - Φ(0,h))
	 * where A(0,h) is the Fresnel term (Equation 12).
	 */
	Float NrIntegrand(Float beta, Float wiDotWo, Float phi, Float h) const;

	/**
	 * @brief Computes the integrand for the azimuthal internal lobes (p>=1), Nₚ.
	 * @note For reference/precomputation only. Not used at runtime.
	 *
	 * This function represents the term inside the integral of Equation (10) for p>=1:
	 *   A(p,h) * Dₚ(φ - Φ(p,h))
	 * where A(p,h) is the full attenuation term including absorption (Equation 13).
	 */
	Vector3 NpIntegrand(Float beta, Float cosThetaD, Float phi, int p, Float h) const;

	void precomputeAzimuthalDistributions();
	void precomputeDS();

	void calculateDrivedProperty();

	/**
	 * @brief Generates a pre-integrated lookup table for the dEonHairBSDF.
	 *
	 * This function iterates through incident angles theta_i, calculates the integral of the BSDF
	 * over the entire sphere of outgoing directions, and saves the results to a binary file.
	 * Due to azimuthal symmetry, the table only needs to be indexed by theta_i.
	 *
	 * @param hairBSDF A pointer to the dEonHairBSDF object.
	 * @param filename The name of the output binary file.
	 * @param theta_i_steps The number of steps for the incident angle theta_i (defines table size).
	 * @param integration_theta_o_steps The number of sampling steps for theta_o for the numerical integration.
	 * @param integration_phi_o_steps The number of sampling steps for phi_o for the numerical integration.
	 */
	void generateIntegratedHairBSDFTable(
		const dEonHairBSDF *hairBSDF,
		const std::string &filename,
		int theta_i_steps = 90,				// 90 points for theta_i from -90 to +90 deg
		int integration_theta_o_steps = 90, // Integration resolution
		int integration_phi_o_steps = 180	// Integration resolution
	);

	/**
	 * @brief Performs numerical integration of the BSDF over the ENTIRE SPHERE for a given incident direction.
	 * This is crucial for transmissive materials like hair, which scatter light in all directions.
	 *
	 * @param hairBSDF The BSDF model to evaluate.
	 * @param wi The incident direction vector.
	 * @param integration_theta_steps The number of steps for theta_o during integration.
	 * @param integration_phi_steps The number of steps for phi_o during integration.
	 * @return The integrated Spectrum value (Directional-Hemispherical Reflectance).
	 */
	Spectrum integrateBSDF(const dEonHairBSDF *hairBSDF, const Vector3f &wi,
						   int integration_theta_steps, int integration_phi_steps);

	Spectrum sampleWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, const Point2 &sample, Float &eta, BSDF::EBSDFType &sampledType, Sampler *sampler) const;
	Float pdfWithLocalWiWo(const Vector3 &wi, const Vector3 &wo) const;
	Spectrum evalWithLocalWiWo(const Vector3 &wi, const Vector3 &wo) const;
};

MTS_NAMESPACE_END

#endif