#ifndef __HAIR_VOLUME_H
#define __HAIR_VOLUME_H

#include <memory>
#include <mitsuba/core/distrib.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/scene.h>
#include <src/phase/dEonPhase.h>
#include <vector>

MTS_NAMESPACE_BEGIN

/**
 * Hair volume integrator.
 *
 * Traces surface paths through explicit hair geometry and evaluates an
 * approximating participating medium after a path leaves a hair curve. The
 * medium event contributes direct illumination using the octave-based
 * multiple-scattering approximation.
 */
class HairVolumeIntegrator : public MonteCarloIntegrator {
  public:
	HairVolumeIntegrator(const Properties &props);
	HairVolumeIntegrator(Stream *stream, InstanceManager *manager);

	void serialize(Stream *stream, InstanceManager *manager) const override;
	void addChild(const std::string &name, ConfigurableObject *child) override;

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const override;

	std::string toString() const override;

	MTS_DECLARE_CLASS()

  protected:
	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	/// Build the sampling distribution used for stochastic octave selection.
	void configureOctaves();

	/// March through the hair medium and estimate segment transmittance.
	Spectrum evalTransmittanceManual(const Point &pStart, const Point &pEnd,
									 const Medium *medium, Float maxDist, Float time, Float extinctionScale) const;

	RayMaskType getRayMask(int depth) const;

	/// Evaluate the contribution scale for one multiple-scattering octave.
	Spectrum evalOctaveWeights(
		int octaveIndex,
		Float scaleG,
		const Spectrum &localAlbedo,
		const Vector3 &rawPhaseRGB,
		Float isotropicVal
	) const;

	/// Sample direct illumination at a medium interaction.
	Spectrum sampleMediumDirectLight(RadianceQueryRecord &rRec, const RayDifferential &ray,
									 const MediumSamplingRecord &mRec, const Spectrum &pathThroughput) const;

  protected:
	Float m_stepSize;
	bool m_jitter;

	int m_msOctaves;
	Float m_msA, m_msB, m_msC;
	Float m_gamma;
	std::unique_ptr<DistribTable> m_octaveDist;

	/// If true, sum every octave. Otherwise, sample one octave and divide by its PDF.
	bool m_deterministicOctaves;

	/// Profiling/debug switches retained from the research implementation.
	bool m_profileDisableWoodcock = false;
	bool m_profileForceScatterAtSurface = false;
	bool m_profileDisableLightMarch = false;
	bool m_profileDisableDeonPhaseEval = false;
	bool m_profileClampToNextHit = false;

	Medium *m_hairMedium = nullptr;
};

MTS_NAMESPACE_END

#endif /* __HAIR_VOLUME_H */
