#include "hairvolume.h"
#include <algorithm>
#include <cmath>

#include <mitsuba/core/profiler.h>
#include <mitsuba/render/renderjob.h>

MTS_NAMESPACE_BEGIN

HairVolumeIntegrator::HairVolumeIntegrator(const Properties &props)
	: MonteCarloIntegrator(props) {
	m_stepSize = props.getFloat("stepSize", 0.1f);
	m_jitter = props.getBoolean("jitter", true);

	m_msOctaves = props.getInteger("octaves", 1);
	m_msA = props.getFloat("a", 0.5f);
	m_msB = props.getFloat("b", 0.5f);
	m_msC = props.getFloat("c", 0.5f);
	m_gamma = props.getFloat("gamma", 1.f);

	m_deterministicOctaves = props.getBoolean(
		props.hasProperty("deterministicOctaves") ? "deterministicOctaves" : "deterministicOz", true);
	m_profileDisableWoodcock = props.getBoolean("profileDisableWoodcock", false);
	m_profileForceScatterAtSurface = props.getBoolean("profileForceScatterAtSurface", false);
	m_profileDisableLightMarch = props.getBoolean("profileDisableLightMarch", false);
	m_profileDisableDeonPhaseEval = props.getBoolean("profileDisableDeonPhaseEval", false);
	m_profileClampToNextHit = props.getBoolean("profileClampToNextHit", false);

	configureOctaves();
}

HairVolumeIntegrator::HairVolumeIntegrator(Stream *stream, InstanceManager *manager)
	: MonteCarloIntegrator(stream, manager) {
	m_stepSize = stream->readFloat();
	m_jitter = stream->readBool();
	m_msOctaves = stream->readInt();
	m_msA = stream->readFloat();
	m_msB = stream->readFloat();
	m_msC = stream->readFloat();
	m_gamma = stream->readFloat();
	m_deterministicOctaves = stream->readBool();
	m_profileDisableWoodcock = stream->readBool();
	m_profileForceScatterAtSurface = stream->readBool();
	m_profileDisableLightMarch = stream->readBool();
	m_profileDisableDeonPhaseEval = stream->readBool();
	m_profileClampToNextHit = stream->readBool();

	configureOctaves();
}

void HairVolumeIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
	MonteCarloIntegrator::serialize(stream, manager);
	stream->writeFloat(m_stepSize);
	stream->writeBool(m_jitter);
	stream->writeInt(m_msOctaves);
	stream->writeFloat(m_msA);
	stream->writeFloat(m_msB);
	stream->writeFloat(m_msC);
	stream->writeFloat(m_gamma);
	stream->writeBool(m_deterministicOctaves);
	stream->writeBool(m_profileDisableWoodcock);
	stream->writeBool(m_profileForceScatterAtSurface);
	stream->writeBool(m_profileDisableLightMarch);
	stream->writeBool(m_profileDisableDeonPhaseEval);
	stream->writeBool(m_profileClampToNextHit);
}

void HairVolumeIntegrator::addChild(const std::string &name, ConfigurableObject *child) {
	if (child->getClass()->derivesFrom(MTS_CLASS(Medium))) {
		Medium *medium = static_cast<Medium *>(child);
		if (m_hairMedium == nullptr)
			m_hairMedium = medium;
	} else {
		MonteCarloIntegrator::addChild(name, child);
	}
}

void HairVolumeIntegrator::configureOctaves() {
	if (m_msOctaves < 1)
		Log(EError, "The 'octaves' parameter must be at least 1");

	std::vector<Float> weights(m_msOctaves);
	for (int i = 0; i < m_msOctaves; ++i)
		weights[i] = std::pow(m_msB, (Float)i);

	m_octaveDist = std::make_unique<DistribTable>(&weights[0], m_msOctaves);
}

Spectrum HairVolumeIntegrator::evalTransmittanceManual(const Point &pStart, const Point &pEnd,
													   const Medium *medium, Float maxDist, Float time, Float extinctionScale) const {
	if (!medium || maxDist < Epsilon)
		return Spectrum(1.0f);

	Vector d = pEnd - pStart;
	Float len = d.length();
	Vector dNorm = d / len;

	// Avoid marching slightly past the sampled light position.
	Float targetDist = std::min(len, maxDist * (1 - ShadowEpsilon));

	Spectrum opticalDepth(0.0f);
	Float t = 0.0f;

	Float shadowStep = m_stepSize;

	while (t < targetDist) {
		Float step = std::min(shadowStep, targetDist - t);

		// Midpoint integration is sufficient for the short shadow segments
		// evaluated by this approximation.
		Point p = pStart + dNorm * (t + step * 0.5f);

		Float sT;
		Vector3 orientation;
		medium->evalDensityWithDirection(p, dNorm, sT);

		if (opticalDepth.max() > 10.0f) {
			return Spectrum(0.0f);
		}

		opticalDepth += Spectrum(sT) * step * extinctionScale;
		t += step;
	}

	return (-opticalDepth).exp();
}

RayMaskType HairVolumeIntegrator::getRayMask(int /*depth*/) const {
	return ERayMaskCommon::EDefault | EHairRayMask::EHairGeometry;
}

// Evaluate the albedo and phase scaling for one multiple-scattering octave.
// The caller precomputes the expensive local albedo and raw phase terms once
// per medium sample, then reuses them across octaves.
Spectrum HairVolumeIntegrator::evalOctaveWeights(
	int octaveIndex,
	Float scaleG,
	const Spectrum &localAlbedo,
	const Vector3 &rawPhaseRGB,
	Float isotropicVal
) const {
	Spectrum contribution;
	Spectrum effectiveAlbedo = localAlbedo * m_gamma;
	for (int k = 0; k < SPECTRUM_SAMPLES; ++k) {
		contribution[k] = std::pow(std::clamp(effectiveAlbedo[k], 0.f, 0.9999f), (Float)octaveIndex);
	}

	if (contribution.isZero())
		return Spectrum(0.f);

	Float interpT = 1.0f - scaleG;

	// Higher octaves are blended toward isotropic scattering by reducing the
	// effective phase eccentricity.
	Vector3 isoVec(isotropicVal);
	Vector3 mixedPhaseVec = math::lerp(interpT, rawPhaseRGB, isoVec);

	Spectrum phaseVal;
	phaseVal.fromLinearRGB(mixedPhaseVec.x, mixedPhaseVec.y, mixedPhaseVec.z);

	return contribution * phaseVal;
}

Spectrum HairVolumeIntegrator::sampleMediumDirectLight(
	RadianceQueryRecord &rRec, const RayDifferential &ray,
	const MediumSamplingRecord &mRec, const Spectrum &pathThroughput) const {

	ManualTimer ssTimer(ProfileCategory::VolumetricSingleScatterTransmittance);
	ssTimer.start();

	Spectrum Li_direct(0.f);

	DirectSamplingRecord dRec(mRec.p, mRec.time);
	Spectrum emitterVal;
	{
		ScopeTimer lightSampTimer(ProfileCategory::VolumetricLightSampling);
		emitterVal = rRec.scene->sampleEmitterDirect(dRec, rRec.nextSample2D(), false);

		if (emitterVal.isZero()) {
			ssTimer.stop();
			return Li_direct;
		}

	}

	// Precompute the medium terms shared by all octave evaluations.
	Spectrum localAlbedo = (mRec.sigmaA + mRec.sigmaS).isZero()
							   ? Spectrum(0.f)
							   : (mRec.sigmaS / (mRec.sigmaA + mRec.sigmaS));

	Float sample1D = rRec.nextSample1D();
	PhaseFunctionSamplingRecord pRec(mRec, -ray.d, dRec.d);

	const PhaseFunction *phaseFunction = mRec.medium->getPhaseFunction();
	const dEonPhaseFunction *deonPhase = dynamic_cast<const dEonPhaseFunction *>(phaseFunction);

	const SpectrumPhaseFunction *spectrumPhase = nullptr;
	if (phaseFunction->getClass()->derivesFrom(MTS_CLASS(SpectrumPhaseFunction))) {
		spectrumPhase = static_cast<const SpectrumPhaseFunction *>(phaseFunction);
	}

	Float isotropicVal = warp::squareToUniformSpherePdf();

	Vector3 rawPhaseRGB;
	if (deonPhase) {
		// Keep directional d'Eon albedo even when phase evaluation is disabled.
		localAlbedo = deonPhase->albedoDir(ray.d, mRec.orientation);
	}

	if (m_profileDisableDeonPhaseEval) {
		rawPhaseRGB = Vector3(isotropicVal);
	} else if (deonPhase) {
		rawPhaseRGB = deonPhase->evalRGB(pRec);
	} else if (spectrumPhase) {
		rawPhaseRGB = spectrumPhase->evalRGB(pRec, sample1D);
	} else {
		Float phaseFloat = phaseFunction->eval(pRec);
		rawPhaseRGB = Vector3(phaseFloat);
	}

	Spectrum radianceCommon = pathThroughput * emitterVal;

	if (m_deterministicOctaves) {
		// Base transmittance on scatter->light segment.
		Spectrum baseVolTrans = m_profileDisableLightMarch
									? Spectrum(1.0f)
									: evalTransmittanceManual(mRec.p, dRec.p, rRec.medium, dRec.dist, ray.time, 1.f);

		ssTimer.stop();
		{
			ScopeTimer msScaleTimer(ProfileCategory::VolumetricMultiScatterScaleLoop);
			for (int i = 0; i < m_msOctaves; ++i) {
				Float scaleT = std::pow(m_msA, (Float)i);
				Float scaleG = std::pow(m_msC, (Float)i);

				Spectrum volTrans;
				volTrans.fromLinearRGB(
					std::pow(baseVolTrans[0], scaleT),
					std::pow(baseVolTrans[1], scaleT),
					std::pow(baseVolTrans[2], scaleT));

				Spectrum weights = evalOctaveWeights(i, scaleG, localAlbedo, rawPhaseRGB, isotropicVal);

				if (!weights.isZero()) {
					Li_direct += radianceCommon * weights * volTrans;
				}
			}
		}

	} else {
		Float octaveSampleU = rRec.nextSample1D();
		Float octavePdf;
		uint32_t octaveIndex = m_octaveDist->sample(octaveSampleU, octavePdf);


		if (octavePdf > 0) {
			Float scaleT = std::pow(m_msA, (Float)octaveIndex);
			Float scaleG = std::pow(m_msC, (Float)octaveIndex);

			Spectrum volTrans = m_profileDisableLightMarch
									? Spectrum(1.0f)
									: evalTransmittanceManual(mRec.p, dRec.p, rRec.medium, dRec.dist, ray.time, scaleT);

			ssTimer.stop();

			if (!volTrans.isZero()) {
				ScopeTimer msScaleTimer(ProfileCategory::VolumetricMultiScatterScaleLoop);

				Spectrum weights = evalOctaveWeights(octaveIndex, scaleG, localAlbedo, rawPhaseRGB, isotropicVal);

				if (!weights.isZero()) {
					Li_direct += (radianceCommon * weights * volTrans) / octavePdf;
				}
			}
		}
	}

	return Li_direct;
}

Spectrum HairVolumeIntegrator::Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
	ScopeTimer shadingTimer(ProfileCategory::Shading);

	const Scene *scene = rRec.scene;
	Intersection &its = rRec.its;
	RayDifferential ray(r);
	ray.mask = ERayMaskCommon::EDefault | EHairRayMask::EHairGeometry;
	Spectrum Li(0.0f);
	bool scattered = false;

	{
		ScopeTimer intersectionTimer(ProfileCategory::FirstIntersection);
		rRec.rayIntersect(ray);
	}
	ray.mint = Epsilon;

	Spectrum throughput(1.0f);
	Float eta = 1.0f;

	while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
		if (!its.isValid()) {
			ScopeTimer evalEnvironmentTimer(ProfileCategory::EvalEnvironment);
			if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered))
				Li += throughput * scene->evalEnvironment(ray);
			break;
		}

		const BSDF *bsdf = its.getBSDF(ray);

		if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && (!m_hideEmitters || scattered))
			Li += throughput * its.Le(-ray.d);

		if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
			Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

		if ((rRec.depth >= m_maxDepth && m_maxDepth > 0) || (m_strictNormals && dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0)) {
			break;
		}

		// Surface next-event estimation.
		DirectSamplingRecord dRec(its);

		if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
			(bsdf->getType() & BSDF::ESmooth)) {
			ScopeTimer surfaceDirectTimer(ProfileCategory::SurfaceDirectLighting);

			Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D(), ERayMaskCommon::EDefault | EHairRayMask::EHairGeometry, true);
			if (!value.isZero()) {
				const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

				BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

				bRec.sampler = rRec.sampler;
				bRec.scene = scene;
				const Spectrum bsdfVal = bsdf->eval(bRec);

				if (!bsdfVal.isZero() && (!m_strictNormals || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {
					Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
										? bsdf->pdf(bRec)
										: 0;

					Float weight = miWeight(dRec.pdf, bsdfPdf);
					Li += throughput * value * bsdfVal * weight;
				}
			}
		}

		// Continue the surface path using BSDF sampling.
		bool hitEmitter = false;
		BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
		bRec.scene = rRec.scene;
		Spectrum bsdfWeight = Spectrum(0.f);
		Float bsdfPdf;
		Spectrum value;

		bool previousItsIsCurve = its.shape->isCurve();
		{
			ScopeTimer surfaceSampleTiemr(ProfileCategory::SurfaceSampling);

			bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
			if (bsdfWeight.isZero())
				break;

			scattered |= bRec.sampledType != BSDF::ENull;

			const Vector wo = its.toWorld(bRec.wo);
			Float woDotGeoN = dot(its.geoFrame.n, wo);
			if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				break;

			Point nextRayOrigin = its.p;
			if (its.shape->isCurve()) {
				nextRayOrigin = its.shape->getNextRayOrigin(its, wo);
			}

			ray = Ray(nextRayOrigin, wo, ray.time);
			ray.mask = ERayMaskCommon::EDefault | EHairRayMask::EHairGeometry;
			if (scene->rayIntersect(ray, its)) {
				if (its.isEmitter()) {
					value = its.Le(-ray.d);
					dRec.setQuery(ray, its);
					hitEmitter = true;
				}
			} else {
				const Emitter *env = scene->getEnvironmentEmitter();

				if (env) {
					if (m_hideEmitters && !scattered)
						break;

					value = env->evalEnvironment(ray);
					if (!env->fillDirectSamplingRecord(dRec, ray))
						break;
					hitEmitter = true;
				} else {
					break;
				}
			}
		}

		throughput *= bsdfWeight;
		eta *= bRec.eta;

		if (hitEmitter &&
			(rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
			const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ? scene->pdfEmitterDirect(dRec) : 0;
			Li += throughput * value * miWeight(bsdfPdf, lumPdf);
		}

		if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
			break;
		rRec.type = RadianceQueryRecord::ERadianceNoEmission;

		if (rRec.depth++ >= m_rrDepth) {
			Float q = std::min(throughput.max() * eta * eta, (Float)0.95f);
			if (rRec.nextSample1D() >= q)
				break;
			throughput /= q;
		}

		if (previousItsIsCurve) {
			// After leaving a hair curve, sample the approximating volume once
			// and add its direct-lighting contribution.
			ScopeTimer fullVolShading(ProfileCategory::FullVolumetricShading);
			MediumSamplingRecord mRec;
			if (m_hairMedium == nullptr)
				Log(EError, "The hairvolume integrator requires a hair medium child");
			rRec.medium = m_hairMedium;
			Ray mediumRay(ray);
			if (its.isValid() && m_profileClampToNextHit)
				mediumRay.maxt = its.t;

			bool sampleDistanceResult = false;
			if (m_profileForceScatterAtSurface) {
				// Falcor parity: bypass Woodcock and force scatter at the surface point.
				mRec.t = 0.0f;
				mRec.p = mediumRay.o;
				mRec.time = mediumRay.time;
				mRec.transmittance = Spectrum(1.0f);
				mRec.pdfFailure = 1.0f;
				mRec.pdfSuccess = 1.0f;
				mRec.pdfSuccessRev = 1.0f;
				mRec.variance = 0.0f;
				Spectrum sigmaT;
				rRec.medium->evalCoefficientsWithDirection(mRec.p, mediumRay.d, mRec.sigmaS, sigmaT, mRec.orientation);
				mRec.sigmaA = sigmaT - mRec.sigmaS;
				mRec.medium = rRec.medium;
				mRec.profileAcceptRatio = 1.0f;
				mRec.profileMaxRatio = 1.0f;
				sampleDistanceResult = true;
			} else if (!m_profileDisableWoodcock) {
				ScopeTimer deltaTrackingTimer(ProfileCategory::DeltaTracking);
				sampleDistanceResult = rRec.medium->sampleDistance(mediumRay, mRec, rRec.sampler);
			}
			if (sampleDistanceResult) {
				if (!m_profileForceScatterAtSurface) {
					throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;
				}

				bool maxDepthReached = (rRec.depth >= m_maxDepth && m_maxDepth != -1);
				if (maxDepthReached)
					break;

				if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
					Li += sampleMediumDirectLight(rRec, ray, mRec, throughput);
				}
			}
			break;
		}
	}

	return Li;
}

std::string HairVolumeIntegrator::toString() const {
	std::ostringstream oss;
	oss << "HairVolumeIntegrator[" << endl
		<< "  stepSize = " << m_stepSize << "," << endl
		<< "  jitter = " << m_jitter << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(HairVolumeIntegrator, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(HairVolumeIntegrator, "Hair volume integrator");
MTS_NAMESPACE_END
