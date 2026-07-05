#pragma once

#include <mitsuba/core/properties.h>
#include <mitsuba/core/shvector.h>
#include <mitsuba/render/phase.h>
#include <string>
#include <vector>

MTS_NAMESPACE_BEGIN

#ifdef _WIN32 // Windows/MSVC symbol export
#ifdef DEON_PHASE_EXPORTS
#define DEON_API __declspec(dllexport)
#else
#define DEON_API __declspec(dllimport)
#endif
#else // GCC/Clang symbol export
#define DEON_API __attribute__((visibility("default")))
#endif

// Forward declarations to reduce header dependencies
class Sampler;
class Stream;
class InstanceManager;
class PrecomputedAzimuthalLobe;

/**
 * @brief Utility function to write the precomputed integration tables to text files
 *        for debugging and visualization.
 */
void printTablesToFiles(
	const std::vector<Float> &scattering_table,
	const std::vector<Vector3> &full_table,
	const std::string &scattering_filename,
	const std::string &absorption_filename, const std::string &full_filename);

/**
 * @brief Implements the d'Eon model for hair scattering, a physically-based
 *        phase function that simulates the R, TT, and TRT lobes.
 */
class DEON_API dEonPhaseFunction : public PhaseFunction {
  public:
	// Constructors
	dEonPhaseFunction(const Properties &props);
	dEonPhaseFunction(Stream *stream, InstanceManager *manager);
	void addChild(const std::string &name, ConfigurableObject *child) override;
	// Overridden PhaseFunction interface methods
	void configure() override;
	Float eval(const PhaseFunctionSamplingRecord &pRec) const override;
	Vector3 evalRGB(const PhaseFunctionSamplingRecord &pRec) const;
	Float sample(PhaseFunctionSamplingRecord &pRec, Sampler *sampler) const override;
	Float sample(PhaseFunctionSamplingRecord &pRec, Float &pdf, Sampler *sampler) const override;
	Float pdf(const PhaseFunctionSamplingRecord &pRec) const override;
	bool needsDirectionallyVaryingCoefficients() const override;
	Float sigmaDir(Float cosTheta) const override;
	Float sigmaDirMax() const override;
	Float getMeanCosine() const override;
	Spectrum albedoDir(const Vector &d, const Vector &orientation) const override;
	bool needsDirectionallyVaryingAlbedo() const override;

	/**
	 * \brief Precomputes the Hair Scattering Matrix F for a given number of SH bands.
	 *
	 * \param mat       The hair material parameters.
	 * \param bands     The number of SH bands (e.g., 3 for 9 coeffs, 4 for 16 coeffs).
	 * \param nSamples  Number of Monte Carlo samples (e.g., 1000000).
	 * \return          A dynamic matrix of size (bands^2) x (bands^2).
	 */
	DynamicMatrix PrecomputeHairMatrix(int bands, int nSamples) const override;
	MatrixRGB PrecomputeHairMatrixRGB(int bands, int nSamples) const override;

	std::string toString() const override;
	inline std::shared_ptr<std::vector<Vector3>> getTotalScatteredEnergyTable() const {
		if (m_totalScatteredEnergy->size() == 0) {
			Log(EError, "Attempt to get an empty integration table.");
		}
		return m_totalScatteredEnergy;
	}

	inline Properties getProperties() const {
		return m_properties;
	}

	const Vector3 getRGBMask() const {
		return rgbMask;
	}

	const int getMaskCount() const {
		return maskCount;
	}

	void deriveNormalizationFactor(std::shared_ptr<std::vector<Vector3>> totalScatteredEnergy);

	MTS_DECLARE_CLASS()

  protected:
	// Core evaluation and sampling logic using local coordinates
	Float evalWithLocalWiWo(const Vector3 &wi, const Vector3 &wo, Float variance) const;
	Vector3 evalRGBWithLocalWiWo(const Vector3 &wi, const Vector3 &wo, Float variance) const;
	PFPrecalcSample evalBSDF(const Vector3 &wi, const Vector3 &wo) const;
	Float sampleWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const;
	Vector3 sampleRGBWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const;
	PFPrecalcSample sampleWithLocalWiForIntegration(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const;
	Float pdfWithLocalWiWo(const Vector3 &wi, const Vector3 &wo) const;

	// Initialization and precomputation methods
	void calculateDrivedProperty();
	void precomputeAzimuthalDistributions();
	void generateIntegratedHairBSDFTable();

	// Core d'Eon model mathematical functions
	Float sampleM(Float v, Float sinThetaI, Float cosThetaI, Float xi1, Float xi2) const;
	Float g(Float beta, Float theta) const;
	Float D(Float beta, Float phi) const;
	Float Phi(Float gammaI, Float gammaT, int p) const;
	Float I0(Float x) const;
	Float logI0(Float x) const;
	Float M(Float v, Float sinThetaI, Float sinThetaO, Float cosThetaI, Float cosThetaO) const;
	Float fresnel_dielectric(Float cos_theta_i, Float eta);

	// Helper and integration methods
	Vector3f anglesToVector(Float theta, Float phi) const;
	PFPrecalcSample integrateBSDF(const Vector3f &wi) const;
	Float integrateBSDFTest(const Vector3f &wi) const;
	Float lookupNormalizationFactor(Float thetaI) const;
	Vector3 lookupTotalScatteredEnergy(Float thetaI) const;
	Vector3 lookupTotalTransmissiveEnergy(Float thetaI) const;

	// Serialization
	// Helper to generate a unique filename hash based on current parameters
	std::string generateCacheFilename() const;
	// Try to load the table from disk. Returns true if successful.
	bool loadPrecomputedTable(const fs::path &path);
	// Save the current table to disk.
	void savePrecomputedTable(const fs::path &path) const;

  protected:
	/// @brief Index of refraction of the hair fiber. Typically ~1.55 for human hair.
	/// Corresponds to the symbol η (eta) in the paper.
	const Float Eta = 1.55;

	// User-configurable appearance parameters
	/// @brief Longitudinal shift of scattering lobes due to cuticle angle, in degrees.
	Float _scaleAngleDeg;

	/// @brief Ratio of eumelanin (brown/black) to pheomelanin (red/yellow). Range: [0, 1].
	/// Used to procedurally calculate the absorption coefficient _sigmaA.
	Float _melaninRatio;

	/// @brief Overall concentration of melanin pigments. Controls hair darkness.
	/// Also used to procedurally calculate _sigmaA.
	Float _melaninConcentration;

	/// @brief The volume absorption coefficient [m^-1], defining the hair color.
	/// This corresponds to mu_a in the paper and is crucial for calculating
	/// the attenuation term A(p,h) for internal lobes (p >= 1).
	Vector3 _sigmaA;

	/// @brief Isotropic roughness of the hair fiber. Controls the width of highlights.
	/// This value is used to derive per-lobe roughness values (β_p).
	Float _roughness;

	// Internal/derived model parameters
	// These are precomputed in the constructor for efficiency.

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

	/// Parameters for integration and table generation
	int m_theta_i_steps = 480.f;
	int m_integration_theta_steps = 180.f;
	int m_integration_phi_steps = 180.f;
	int m_intergration_sample_cnt = 65536;

	/// Storage for the precomputed table
	std::vector<Float> m_normalizationFactors;
	std::shared_ptr<std::vector<Vector3>> m_totalScatteredEnergy;
	std::shared_ptr<std::vector<Vector3>> m_totalTransmissiveEnergy;
	Properties m_propterties; // Cached properties used for cloning

	Vector3 rgbMask;
	Float maskCount;
	bool maskPhase;
	bool maskAlbedo;
	bool useUniformSphereSampling;
	bool m_overrideSigmaA = false;
	bool m_useTransmissive = true;
};

MTS_NAMESPACE_END
