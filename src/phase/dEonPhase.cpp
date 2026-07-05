#include "dEonPhase.h"
#include <fstream>
#include <functional>
#include <iomanip>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <src/bsdfs/precomputed_marschner.h>

MTS_NAMESPACE_BEGIN

dEonPhaseFunction::dEonPhaseFunction(const Properties &props)
	: PhaseFunction(props) {
	m_totalScatteredEnergy = std::make_shared<std::vector<Vector3>>();
	m_totalTransmissiveEnergy = std::make_shared<std::vector<Vector3>>();

	_roughness = math::clamp(props.getFloat("roughness", 0.5), 0.f, 1.f);
	_melaninConcentration = props.getFloat("melaninConcentration", 0.5);
	_melaninRatio = props.getFloat("melaninRatio", 0.5);
	_scaleAngleDeg = props.getFloat("scaleAngleDeg", 0.f);
	rgbMask = props.getVector("rgbMask", Vector3(1.f));

	maskCount = (std::abs(rgbMask.x) > Epsilon) + (std::abs(rgbMask.y) > Epsilon) + (std::abs(rgbMask.z) > Epsilon);

	maskPhase = props.getBoolean("maskPhase", true);
	maskAlbedo = props.getBoolean("maskAlbedo", true);

	useUniformSphereSampling = props.getBoolean("useUniformSphereSampling", true);

	m_overrideSigmaA = props.getBoolean("overrideSigmaA", false);
	if (m_overrideSigmaA) {
		_sigmaA = props.getVector("sigmaA"); // no default value
	}

	calculateDrivedProperty();
	precomputeAzimuthalDistributions();

	bool calculateIntegration = props.getBoolean("calculateIntegration", true);
	if (calculateIntegration) {
		// 1. Generate a cache filename from the current parameter set.
		std::string cacheName = generateCacheFilename();
		fs::path cacheFilename(cacheName);

		// 2. Resolve the cache path through Mitsuba's file resolver.
		// If it is not found in search paths, resolve() returns the input path.
		const FileResolver *resolver = Thread::getThread()->getFileResolver();
		fs::path resolvedPath = resolver->resolve(cacheFilename);

		bool loaded = false;

		// 3. Load existing cache data when available.
		if (fs::exists(resolvedPath)) {
			loaded = loadPrecomputedTable(resolvedPath);
		}
		if (loaded) {
			deriveNormalizationFactor(m_totalScatteredEnergy);
		}

		// 4. If loading fails (or no file exists), regenerate and save the cache.
		if (!loaded) {
			// Regenerate integration tables.
			generateIntegratedHairBSDFTable();

			// Save to the current working directory (or where resolvedPath points if absolute)
			// Note: If resolve() failed to find it, it usually returns a relative path.
			// We write to the local directory to populate the cache.
			savePrecomputedTable(cacheFilename);
		}

	}
	m_propterties = Properties(props);
	m_useTransmissive = props.getBoolean("useTransmissive", true);
}

void dEonPhaseFunction::addChild(const std::string &name, ConfigurableObject *child) {
	dEonPhaseFunction *otherDEonPhase = dynamic_cast<dEonPhaseFunction *>(child);
	if (otherDEonPhase == nullptr) {
		Log(EError, "Try to add a non dEonPhaseFunction to a dEonPhaseFunction!");
		return;
	}
	deriveNormalizationFactor(otherDEonPhase->getTotalScatteredEnergyTable());
}

dEonPhaseFunction::dEonPhaseFunction(Stream *stream, InstanceManager *manager)
	: PhaseFunction(stream, manager) {
	configure();
}

void dEonPhaseFunction::configure() {
	PhaseFunction::configure();
	m_type = EAnisotropic;
}

Float dEonPhaseFunction::eval(const PhaseFunctionSamplingRecord &pRec) const {
	if (pRec.mRec.orientation.isZero()) {
		return 0.0f;
	}
	Frame frame = Frame::fromS(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);
	Vector wo = frame.toLocal(pRec.wo);
	return evalWithLocalWiWo(wi, wo, pRec.mRec.variance);
}

Float dEonPhaseFunction::evalWithLocalWiWo(const Vector3 &wi, const Vector3 &wo, Float variance) const {
	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);

	Float sinThetaO = wo.x;
	Float cosThetaO = math::trigInverse(sinThetaO);

	Float thetaI = math::safe_asin(sinThetaI);
	Float thetaO = math::safe_asin(sinThetaO);
	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phiI = std::atan2(wi.z, wi.y);
	Float phiO = std::atan2(wo.z, wo.y);
	Float phi = phiO - phiI;
	if (phi < 0.0f)
		phi += 2 * M_PI;

	// Lobe shifts caused by cuticle tilt.
	// Uses the convention from "Importance Sampling for Physically-Based Hair Fiber Models".
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Evaluate longitudinal scattering functions
	Float geoVariance = variance;
	Float MR = M(_vR + geoVariance, std::sin(thetaIR), sinThetaO, std::cos(thetaIR), cosThetaO);
	Float MTT = M(_vTT + geoVariance, std::sin(thetaITT), sinThetaO, std::cos(thetaITT), cosThetaO);
	Float MTRT = M(_vTRT + geoVariance, std::sin(thetaITRT), sinThetaO, std::cos(thetaITRT), cosThetaO);

	Vector3 R = MR * _NR->eval(phi, cosThetaD);
	Vector3 TT = MTT * _NTT->eval(phi, cosThetaD);
	Vector3 TRT = MTRT * _NTRT->eval(phi, cosThetaD);
	Vector3 combinedLobes = R + TT + TRT;

	Vector3 S = (combinedLobes);
	Float C = lookupNormalizationFactor(thetaI);

	Float result = 0.f;
	if (maskPhase) {
		result = CwiseMul(S, rgbMask).sum() / maskCount;
	} else {
		result = S.avg();
	}

	return result / C;
}

Vector3 dEonPhaseFunction::evalRGBWithLocalWiWo(const Vector3 &wi, const Vector3 &wo, Float variance) const {
	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);

	Float sinThetaO = wo.x;
	Float cosThetaO = math::trigInverse(sinThetaO);

	Float thetaI = math::safe_asin(sinThetaI);
	Float thetaO = math::safe_asin(sinThetaO);
	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phiI = std::atan2(wi.z, wi.y);
	Float phiO = std::atan2(wo.z, wo.y);
	Float phi = phiO - phiI;
	if (phi < 0.0f)
		phi += 2 * M_PI;

	// Lobe shifts caused by cuticle tilt.
	// Uses the convention from "Importance Sampling for Physically-Based Hair Fiber Models".
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Evaluate longitudinal scattering functions
	Float geoVariance = variance;
	Float MR = M(_vR + geoVariance, std::sin(thetaIR), sinThetaO, std::cos(thetaIR), cosThetaO);
	Float MTT = M(_vTT + geoVariance, std::sin(thetaITT), sinThetaO, std::cos(thetaITT), cosThetaO);
	Float MTRT = M(_vTRT + geoVariance, std::sin(thetaITRT), sinThetaO, std::cos(thetaITRT), cosThetaO);

	Vector3 R = MR * _NR->eval(phi, cosThetaD);
	Vector3 TT = MTT * _NTT->eval(phi, cosThetaD);
	Vector3 TRT = MTRT * _NTRT->eval(phi, cosThetaD);
	Vector3 combinedLobes = R + TT + TRT;

	Vector3 S = (combinedLobes);
	Vector3 C = lookupTotalScatteredEnergy(thetaI);

	Vector3 result(0.f);
	result = CwiseDiv(S, C);

	return result;
}

Vector3 dEonPhaseFunction::evalRGB(const PhaseFunctionSamplingRecord &pRec) const {
	if (pRec.mRec.orientation.isZero()) {
		return Vector3(0.0f);
	}
	Frame frame = Frame::fromS(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);
	Vector wo = frame.toLocal(pRec.wo);

	return evalRGBWithLocalWiWo(wi, wo, pRec.mRec.variance);
}

PFPrecalcSample dEonPhaseFunction::evalBSDF(const Vector3 &wi, const Vector3 &wo) const {
	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);

	Float sinThetaO = wo.x;
	Float cosThetaO = math::trigInverse(sinThetaO);

	Float thetaI = math::safe_asin(sinThetaI);
	Float thetaO = math::safe_asin(sinThetaO);
	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phiI = std::atan2(wi.z, wi.y);
	Float phiO = std::atan2(wo.z, wo.y);
	Float phi = phiO - phiI;
	if (phi < 0.0f)
		phi += 2 * M_PI;

	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Evaluate longitudinal scattering functions
	Float MR = M(_vR, std::sin(thetaIR), sinThetaO, std::cos(thetaIR), cosThetaO);
	Float MTT = M(_vTT, std::sin(thetaITT), sinThetaO, std::cos(thetaITT), cosThetaO);
	Float MTRT = M(_vTRT, std::sin(thetaITRT), sinThetaO, std::cos(thetaITRT), cosThetaO);

	Vector3 R = MR * _NR->eval(phi, cosThetaD);
	Vector3 TT = MTT * _NTT->eval(phi, cosThetaD);
	Vector3 TRT = MTRT * _NTRT->eval(phi, cosThetaD);
	Vector3 combinedLobes = R + TT + TRT;

	Float scattering = 0.f;
	Vector3 full(0.f);

	Vector3 transmissive = TT + TRT;

	Vector3 S = (combinedLobes);
	if (maskPhase) {
		scattering = CwiseMul(S, rgbMask).sum() / maskCount;
	} else {
		scattering = S.avg();
	}
	full = S;

	return {scattering, full, transmissive};
}

Float dEonPhaseFunction::sample(PhaseFunctionSamplingRecord &pRec,
								Sampler *sampler) const {
	Float pdf = 0.f;
	Float value = this->sample(pRec, pdf, sampler);
	return value;
}

Float dEonPhaseFunction::sample(PhaseFunctionSamplingRecord &pRec,
								Float &pdf, Sampler *sampler) const {
	if (pRec.mRec.orientation.isZero()) {
		return 0.0f;
	}
	Frame frame = Frame::fromS(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);

	Vector wo;
	Float sampleValue = sampleWithLocalWi(wi, wo, pdf, sampler);
	pRec.wo = frame.toWorld(wo);
	return sampleValue;
}

Float dEonPhaseFunction::sampleWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const {
	if (useUniformSphereSampling) {
		wo = warp::squareToUniformSphere(sampler->next2D());
		Float value = this->evalWithLocalWiWo(wi, wo, 0.f);
		pdf = warp::squareToUniformSpherePdf();
		return value / pdf;
	}

	Point2 xiN = sampler->next2D();
	Point2 xiM = sampler->next2D();

	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = math::safe_asin(sinThetaI);

	Float phiI = std::atan2(wi.z, wi.y);

	// Lobe shifts caused by cuticle tilt.
	// Uses the convention from "Importance Sampling for Physically-Based Hair Fiber Models".
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Select one lobe according to its azimuthal importance weight.
	Float weightR = _NR->weight(cosThetaI);
	Float weightTT = _NTT->weight(cosThetaI);
	Float weightTRT = _NTRT->weight(cosThetaI);

	const PrecomputedAzimuthalLobe *lobe;
	Float v;
	Float theta;

	Float target = xiN.x * (weightR + weightTT + weightTRT);
	if (target < weightR) {
		v = _vR;
		theta = thetaIR;
		lobe = _NR.get();
	} else if (target < weightR + weightTT) {
		v = _vTT;
		theta = thetaITT;
		lobe = _NTT.get();
	} else {
		v = _vTRT;
		theta = thetaITRT;
		lobe = _NTRT.get();
	}

	Float sinThetaO = sampleM(v, std::sin(theta), std::cos(theta), xiM.x, xiM.y);
	Float cosThetaO = math::trigInverse(sinThetaO);
	Float thetaO = math::safe_asin(sinThetaO);

	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phi, phiPdf;
	lobe->sample(cosThetaD, xiN.y, phi, phiPdf);
	if (!std::isfinite(phi)) {
		Log(EWarn, "dEonPhaseFunction::sampleWithLocalWi(): non-finite azimuth; returning zero contribution.");
		pdf = 0.f;
		return 0.f;
	}
	Float sinPhiO = std::sin(phi + phiI);
	Float cosPhiO = std::cos(phi + phiI);

	wo = Vector3(
		sinThetaO,
		cosPhiO * cosThetaO,
		sinPhiO * cosThetaO);

	pdf = this->pdfWithLocalWiWo(wi, wo);

	if (pdf > 0.f) {
		return this->evalWithLocalWiWo(wi, wo, 0.f) / pdf;
	} else {
		return 0.f;
	}
}

Vector3 dEonPhaseFunction::sampleRGBWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const {
	if (useUniformSphereSampling) {
		wo = warp::squareToUniformSphere(sampler->next2D());
		Vector3 value = this->evalRGBWithLocalWiWo(wi, wo, 0.f);
		pdf = warp::squareToUniformSpherePdf();
		return value / pdf;
	}

	Point2 xiN = sampler->next2D();
	Point2 xiM = sampler->next2D();

	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = math::safe_asin(sinThetaI);

	Float phiI = std::atan2(wi.z, wi.y);

	// Lobe shifts caused by cuticle tilt.
	// Uses the convention from "Importance Sampling for Physically-Based Hair Fiber Models".
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Select one lobe according to its azimuthal importance weight.
	Float weightR = _NR->weight(cosThetaI);
	Float weightTT = _NTT->weight(cosThetaI);
	Float weightTRT = _NTRT->weight(cosThetaI);

	const PrecomputedAzimuthalLobe *lobe;
	Float v;
	Float theta;

	Float target = xiN.x * (weightR + weightTT + weightTRT);
	if (target < weightR) {
		v = _vR;
		theta = thetaIR;
		lobe = _NR.get();
	} else if (target < weightR + weightTT) {
		v = _vTT;
		theta = thetaITT;
		lobe = _NTT.get();
	} else {
		v = _vTRT;
		theta = thetaITRT;
		lobe = _NTRT.get();
	}

	Float sinThetaO = sampleM(v, std::sin(theta), std::cos(theta), xiM.x, xiM.y);
	Float cosThetaO = math::trigInverse(sinThetaO);
	Float thetaO = math::safe_asin(sinThetaO);

	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phi, phiPdf;
	lobe->sample(cosThetaD, xiN.y, phi, phiPdf);
	if (!std::isfinite(phi)) {
		Log(EWarn, "dEonPhaseFunction::sampleRGBWithLocalWi(): non-finite azimuth; returning zero contribution.");
		pdf = 0.f;
		return Vector3(0.f);
	}
	Float sinPhiO = std::sin(phi + phiI);
	Float cosPhiO = std::cos(phi + phiI);

	wo = Vector3(
		sinThetaO,
		cosPhiO * cosThetaO,
		sinPhiO * cosThetaO);

	pdf = this->pdfWithLocalWiWo(wi, wo);

	if (pdf > 0.f) {
		return this->evalRGBWithLocalWiWo(wi, wo, 0.f) / pdf;
	} else {
		return Vector3(0.f);
	}
}

PFPrecalcSample dEonPhaseFunction::sampleWithLocalWiForIntegration(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const {
	Point2 xiN = sampler->next2D();
	Point2 xiM = sampler->next2D();

	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = math::safe_asin(sinThetaI);

	Float phiI = std::atan2(wi.z, wi.y);

	// Lobe shifts caused by cuticle tilt.
	// Uses the convention from "Importance Sampling for Physically-Based Hair Fiber Models".
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Select one lobe according to its azimuthal importance weight.
	Float weightR = _NR->weight(cosThetaI);
	Float weightTT = _NTT->weight(cosThetaI);
	Float weightTRT = _NTRT->weight(cosThetaI);

	const PrecomputedAzimuthalLobe *lobe;
	Float v;
	Float theta;

	Float target = xiN.x * (weightR + weightTT + weightTRT);
	if (target < weightR) {
		v = _vR;
		theta = thetaIR;
		lobe = _NR.get();
	} else if (target < weightR + weightTT) {
		v = _vTT;
		theta = thetaITT;
		lobe = _NTT.get();
	} else {
		v = _vTRT;
		theta = thetaITRT;
		lobe = _NTRT.get();
	}

	Float sinThetaO = sampleM(v, std::sin(theta), std::cos(theta), xiM.x, xiM.y);
	Float cosThetaO = math::trigInverse(sinThetaO);
	Float thetaO = math::safe_asin(sinThetaO);

	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phi, phiPdf;
	lobe->sample(cosThetaD, xiN.y, phi, phiPdf);
	if (!std::isfinite(phi)) {
		Log(EWarn, "dEonPhaseFunction::sampleWithLocalWiForIntegration(): non-finite azimuth; returning zero contribution.");
		pdf = 0.f;
		return PFPrecalcSample();
	}
	Float sinPhiO = std::sin(phi + phiI);
	Float cosPhiO = std::cos(phi + phiI);

	wo = Vector3(
		sinThetaO,
		cosPhiO * cosThetaO,
		sinPhiO * cosThetaO);

	pdf = this->pdfWithLocalWiWo(wi, wo);
	if (!std::isfinite(pdf)) {
		Log(EWarn, "dEonPhaseFunction::sampleWithLocalWiForIntegration(): non-finite PDF; returning zero contribution.");
		pdf = 0.f;
		return PFPrecalcSample();
	}

	if (pdf > 0.f) {
		auto result = this->evalBSDF(wi, wo) / pdf;
		if (!std::isfinite(result.m_normalizationFactor)) {
			Log(EWarn, "dEonPhaseFunction::sampleWithLocalWiForIntegration(): non-finite sample value; returning zero contribution.");
			return PFPrecalcSample();
		}
		return result;
	} else {
		return PFPrecalcSample();
	}
}

Float dEonPhaseFunction::pdf(const PhaseFunctionSamplingRecord &pRec) const {
	if (useUniformSphereSampling) {
		return warp::squareToUniformSpherePdf();
	}
	Frame frame = Frame::fromS(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi);
	Vector wo = frame.toLocal(pRec.wo);
	return pdfWithLocalWiWo(wi, wo);
}

Float dEonPhaseFunction::pdfWithLocalWiWo(const Vector3 &wi, const Vector3 &wo) const {
	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float sinThetaO = wo.x;
	Float cosThetaO = math::trigInverse(sinThetaO);

	Float thetaI = math::safe_asin(sinThetaI);
	Float thetaO = math::safe_asin(sinThetaO);
	Float thetaD = (thetaO - thetaI) * 0.5f;
	Float cosThetaD = std::cos(thetaD);

	Float phiI = std::atan2(wi.z, wi.y);
	Float phiO = std::atan2(wo.z, wo.y);
	Float phi = phiO - phiI;
	if (phi < 0.0f)
		phi += 2 * M_PI;

	if (!std::isfinite(phi)) {
		Log(EError, "Phi not finite. Phi: %f, PhiI: %f, PhiO: %f", phi, phiI, phiO);
	}

	// Lobe shifts caused by cuticle tilt.
	// Uses the convention from "Importance Sampling for Physically-Based Hair Fiber Models".
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;
	Float weightR = _NR->weight(cosThetaI);
	Float weightTT = _NTT->weight(cosThetaI);
	Float weightTRT = _NTRT->weight(cosThetaI);
	Float weightSum = weightR + weightTT + weightTRT;

	Float pdfR = weightR * M(_vR, std::sin(thetaIR), sinThetaO, std::cos(thetaIR), cosThetaO);
	Float pdfTT = weightTT * M(_vTT, std::sin(thetaITT), sinThetaO, std::cos(thetaITT), cosThetaO);
	Float pdfTRT = weightTRT * M(_vTRT, std::sin(thetaITRT), sinThetaO, std::cos(thetaITRT), cosThetaO);

	Float result = (1.0f / weightSum) *
				   (pdfR * _NR->pdf(phi, cosThetaD) + pdfTT * _NTT->pdf(phi, cosThetaD) + pdfTRT * _NTRT->pdf(phi, cosThetaD));

	return result;
}

bool dEonPhaseFunction::needsDirectionallyVaryingCoefficients() const {
	return true;
}

Float dEonPhaseFunction::sigmaDir(Float cosTheta) const {
	Float sinTheta = math::trigInverse(cosTheta);
	return sinTheta;
}

Float dEonPhaseFunction::sigmaDirMax() const {
	return 1.f;
}

Float dEonPhaseFunction::getMeanCosine() const {
	return 0;
}

Spectrum dEonPhaseFunction::albedoDir(const Vector &d, const Vector &orientation) const {
	Frame frame = Frame::fromS(orientation);
	Vector wi = frame.toLocal(-d);
	Float sinThetaI = wi.x;
	Float thetaI = math::safe_asin(sinThetaI);
	Vector3 scattered = lookupTotalScatteredEnergy(thetaI);
	Vector3 transmissive = lookupTotalTransmissiveEnergy(thetaI);
	Spectrum result;

	if (m_useTransmissive) {
		result.fromLinearRGB(transmissive.x, transmissive.y, transmissive.z);
	} else {
		result.fromLinearRGB(scattered.x, scattered.y, scattered.z);
	}
	return result;
}
// #pragma optimize("", on)

bool dEonPhaseFunction::needsDirectionallyVaryingAlbedo() const {
	return true;
}

DynamicMatrix dEonPhaseFunction::PrecomputeHairMatrix(int bands, int nSamples) const {
	const int numCoeffs = bands * bands;

	// 1. Initialize Result Matrix (Dynamic Size)
	// We start with Zero
	DynamicMatrix F = DynamicMatrix::Zero(numCoeffs, numCoeffs);

	// 2. Prepare reusable storage for SH coefficients
	// Using Eigen vectors allows us to use vector/matrix math later
	DynamicVector sh_in_vals(numCoeffs);
	DynamicVector sh_out_vals(numCoeffs);

	// Integration domain: Area of Input Sphere * Area of Output Sphere
	const Float domainSize = (4.0 * M_PI);

	auto pluginManager = PluginManager::getInstance();
	Properties samplerProperties("ldsampler");
	samplerProperties.setInteger("sampleCount", nSamples);
	ref<ConfigurableObject> samplerObject = pluginManager->createObject(MTS_CLASS(Sampler), samplerProperties);
	ref<Sampler> sampler = dynamic_cast<Sampler *>(samplerObject.get());
	sampler->generate(Point2i(-1));

	// 3. Monte Carlo Integration Loop
	for (int s = 0; s < nSamples; ++s) {
		// A. Sample Directions (Uniformly on sphere)
		Vector wi = warp::squareToUniformSphere(sampler->next2D());
		Vector wo;
		// B. Evaluate BSDF and Geometry Term
		// f_s(...)
		// Combined scalar weight for this sample
		Float pdf = 0.f;
		Float weight_bsdf = this->sampleWithLocalWi(wi, wo, pdf, sampler);

		// If weight is effectively zero, skip SH evaluation to save time
		if (std::abs(weight_bsdf) < 1e-8f) {
			sampler->advance();
			continue;
		}

		// C. Evaluate SH Basis functions for both directions
		EvalSH_All_Eigen(bands, wi, sh_in_vals);
		EvalSH_All_Eigen(bands, wo, sh_out_vals);

		// D. Accumulate into Matrix F
		// We want to add: weight * (sh_out_vals * sh_in_vals^T)
		// Mathematically: F_jk += weight * Y_j(wo) * Y_k(wi)
		// In Linear Algebra: Outer Product

		// Eigen optimization: noalias() prevents creation of a temporary matrix
		F.noalias() += weight_bsdf * (sh_out_vals * sh_in_vals.transpose());
		sampler->advance();
	}

	// 4. Normalize the integral
	// Integral approx = (Volume / N) * Sum
	Float normalizationFactor = domainSize / static_cast<Float>(nSamples);
	F *= normalizationFactor;

	return F;
}

MatrixRGB dEonPhaseFunction::PrecomputeHairMatrixRGB(int bands, int nSamples) const {
	const int numCoeffs = bands * bands;
	const Float domainSize = (4.0 * M_PI);

	// 1. Initialize Result Matrices
	MatrixRGB result;
	result.r = DynamicMatrix::Zero(numCoeffs, numCoeffs);
	result.g = DynamicMatrix::Zero(numCoeffs, numCoeffs);
	result.b = DynamicMatrix::Zero(numCoeffs, numCoeffs);

	// 2. Prepare reusable storage
	DynamicVector sh_in_vals(numCoeffs);
	DynamicVector sh_out_vals(numCoeffs);

	// Sampler setup
	auto pluginManager = PluginManager::getInstance();
	Properties samplerProperties("ldsampler");
	samplerProperties.setInteger("sampleCount", nSamples);
	ref<ConfigurableObject> samplerObject = pluginManager->createObject(MTS_CLASS(Sampler), samplerProperties);
	ref<Sampler> sampler = dynamic_cast<Sampler *>(samplerObject.get());
	sampler->generate(Point2i(-1));

	// Progress reporting interval (every 10%, clamped for small sample counts).
	int reportInterval = std::max(1, nSamples / 10);

	SLog(EInfo, "Starting Hair Matrix Precomputation (%d samples, %d bands)...", nSamples, bands);

	// 3. Monte Carlo Integration Loop
	for (int s = 0; s < nSamples; ++s) {
		// A. Sample Directions
		Vector wi = warp::squareToUniformSphere(sampler->next2D());
		Vector wo;
		Float pdf = 0.f;

		// B. Evaluate BSDF (RGB Version)
		Vector3 weight_spec = this->sampleRGBWithLocalWi(wi, wo, pdf, sampler);

		// Optimization: Skip if energy is negligible in all channels
		Float max_weight = weight_spec.max();
		if (max_weight < 1e-8f) {
			// Even if we skip, we still check if we need to report progress
			if ((s + 1) % reportInterval == 0) {
				Float progress = 100.0f * static_cast<Float>(s + 1) / static_cast<Float>(nSamples);
				SLog(EInfo, "Precomputing... %.0f%% complete", progress);
			}
			sampler->advance();
			continue;
		}

		Float w_r = weight_spec.x;
		Float w_g = weight_spec.y;
		Float w_b = weight_spec.z;

		// C. Evaluate SH Basis
		EvalSH_All_Eigen(bands, wi, sh_in_vals);
		EvalSH_All_Eigen(bands, wo, sh_out_vals);

		// D. Accumulate into Matrices (Per Channel)

		// Red Channel
		if (w_r > 1e-9f) {
			result.r.noalias() += w_r * (sh_out_vals * sh_in_vals.transpose());
		}

		// Green Channel
		if (w_g > 1e-9f) {
			result.g.noalias() += w_g * (sh_out_vals * sh_in_vals.transpose());
		}

		// Blue Channel
		if (w_b > 1e-9f) {
			result.b.noalias() += w_b * (sh_out_vals * sh_in_vals.transpose());
		}

		// Report progress every reportInterval samples.
		if ((s + 1) % reportInterval == 0) {
			Float progress = 100.0f * static_cast<Float>(s + 1) / static_cast<Float>(nSamples);
			SLog(EInfo, "Precomputing... %.0f%% complete (%d / %d samples)", progress, s + 1, nSamples);
		}

		sampler->advance();
	}

	// 4. Normalize the integrals
	Float normalizationFactor = domainSize / static_cast<Float>(nSamples);

	result.r *= normalizationFactor;
	result.g *= normalizationFactor;
	result.b *= normalizationFactor;

	SLog(EInfo, "Hair Matrix Precomputation Finished.");

	return result;
}

std::string dEonPhaseFunction::toString() const {
	std::ostringstream oss;
	oss << "dEonPhase[" << endl
		<< "]";
	return oss.str();
}

void dEonPhaseFunction::calculateDrivedProperty() {
	// The two melanin parameters allow easy reproduction of physical hair colors
	// based on the mixture of two pigments, eumelanin and pheomelanin, found in human hair
	// These RGB absorption values are taken from "An Energy-Conserving Hair Reflectance Model"
	const Vector3 eumelaninSigmaA = Vector3(0.419f, 0.697f, 1.37f);
	const Vector3 pheomelaninSigmaA = Vector3(0.187f, 0.4f, 1.05f);
	if (!m_overrideSigmaA) {
		_sigmaA = _melaninConcentration * math::lerp(_melaninRatio, eumelaninSigmaA, pheomelaninSigmaA);
	}
	// map roughness[0,1] to betaR[M_PI/36, M_PI/18]([5 degs to 10 degs])
	//_betaR = M_PI / 36.f + M_PI / 36.f * _roughness;
	_betaR = std::max(M_PI * 0.5f * _roughness, 0.04f);
	_betaTT = _betaR * 0.5f;
	_betaTRT = _betaR * 2.0f;

	_vR = _betaR * _betaR;
	_vTT = _betaTT * _betaTT;
	_vTRT = _betaTRT * _betaTRT;
	_scaleAngleRad = _scaleAngleDeg / 180 * M_PI;
}

Float dEonPhaseFunction::sampleM(Float v, Float sinThetaI, Float cosThetaI, Float xi1, Float xi2) const {
	// Sample the polar angle component using the stable inversion formula.
	// This gives us a new `cosTheta` relative to the fiber tangent.
	/*Float cosTheta = 1.0f + v * std::log(xi1 + (1.0f - xi1) * std::exp(-2.0f / v));*/
	Float safeExp = std::exp(-2.0f / std::max(v, 1e-6f));
	Float inner = xi1 + (1.0f - xi1) * safeExp;
	inner = std::max(inner, 1e-12f);

	Float cosTheta = 1.0f + v * std::log(inner);

	Float sinTheta = math::trigInverse(cosTheta);

	// Sample the azimuthal angle component uniformly.
	Float cosPhi = std::cos(2 * M_PI * xi2);

	// Combine the sampled components with the incident angle using spherical trigonometry
	// to get the final outgoing longitudinal angle sin(θ_o).
	if (!std::isfinite(cosTheta) || !std::isfinite(cosTheta) || !std::isfinite(cosPhi)) {
		Log(EWarn, "sampleM(): non-finite intermediate value; returning zero longitudinal angle.");
		return 0.f;
	}

	return -cosTheta * sinThetaI + sinTheta * cosPhi * cosThetaI;
}

Float dEonPhaseFunction::g(Float beta, Float theta) const {
	return std::exp(-theta * theta / (2.0f * beta * beta)) / (math::safe_sqrt(2.0f * M_PI) * beta);
}

Float dEonPhaseFunction::D(Float beta, Float phi) const {
	Float result = 0.0f;
	Float delta;
	Float shift = 0.0f;
	do {
		// Sum Gaussians centered at +shift and -(shift + 2pi), etc.
		// This covers the wrapping nature of the angle phi.
		delta = g(beta, phi + shift) + g(beta, phi - shift - 2 * M_PI);
		result += delta;
		shift += 2 * M_PI;
	} while (delta > 1e-4f); // Stop when terms are too small to contribute
	return result;
}

Float dEonPhaseFunction::Phi(Float gammaI, Float gammaT, int p) const {
	return 2.0f * p * gammaT - 2.0f * gammaI + p * M_PI;
}

Float dEonPhaseFunction::I0(Float x) const {
	Float result = 1.0f;
	Float xSq = x * x;
	Float xi = xSq;
	Float denom = 4.0f; // Represents 4^i * (i!)^2, starting with i=1
	for (int i = 1; i <= 10; ++i) {
		result += xi / denom;
		xi *= xSq;
		// Update denominator for the next term: 4*(i+1)^2
		denom *= 4.0f * Float((i + 1) * (i + 1));
	}
	return result;
}

Float dEonPhaseFunction::logI0(Float x) const {
	if (x > 12.0f) {
		// Asymptotic expansion for log(I₀(x)) for large x:
		// log(I₀(x)) ≈ x - 0.5*log(2πx) + 1/(8x)
		return x + 0.5f * (std::log(1.0f / (2 * M_PI * x)) + 1.0f / (8.0f * x));
	} else {
		return std::log(I0(x));
	}
}

Float dEonPhaseFunction::M(Float v, Float sinThetaI, Float sinThetaO, Float cosThetaI, Float cosThetaO) const {
	// The arguments to the exp and I₀ functions
	Float a = cosThetaI * cosThetaO / v;
	Float b = sinThetaI * sinThetaO / v; // Note: the paper uses sin(-θᵢ) = -sin(θᵢ)

	if (v < 0.1f) {
		// Use a numerically stable evaluation for small roughness 'v'.
		// This combines log-space calculations to prevent underflow/overflow.
		// log(M_p) = -b + log(I₀(a)) - log(2v * sinh(1/v))
		// log(2v*sinh(1/v)) = log(e^(1/v) - e^(-1/v)) ≈ 1/v - log(2) for small v
		// So log(M_p) ≈ -b + log(I₀(a)) - (1/v - log(2))
		return std::exp(-b + logI0(a) - 1.0f / v + 0.6931f + std::log(1.0f / (2.0f * v)));
	} else {
		// Direct implementation of Equation (7), using csch(x) = 1/sinh(x)
		return std::exp(-b) * I0(a) / (2.0f * v * std::sinh(1.0f / v));
	}
}

inline Float dEonPhaseFunction::fresnel_dielectric(Float cos_theta_i, Float eta) {
	Float sin_theta_t_2 = eta * eta * (1.0f - cos_theta_i * cos_theta_i);

	// Total internal reflection
	if (sin_theta_t_2 > 1.0)
		return 1.0;

	Float cos_theta_t = std::sqrt(std::max(1.0f - sin_theta_t_2, 0.0f));

	Float rs = (eta * cos_theta_t - cos_theta_i) / (eta * cos_theta_t + cos_theta_i);
	Float rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

	return 0.5f * (rs * rs + rp * rp);
}

void dEonPhaseFunction::precomputeAzimuthalDistributions() {
	const int Resolution = PrecomputedAzimuthalLobe::AzimuthalResolution;
	std::vector<Vector3> valuesR(Resolution * Resolution);
	std::vector<Vector3> valuesTT(Resolution * Resolution);
	std::vector<Vector3> valuesTRT(Resolution * Resolution);

	const int numPoints = 140;
	std::vector<Float> points(numPoints);
	std::vector<Float> weights(numPoints);
	gaussLegendre(numPoints, points.data(), weights.data());

	// We are precomputing the integral of N_p(phi): Equation (10) by integrating over the fiber's cross-section.
	// The numerical integration method (Gauss-Legendre) gives us sample points 'h' in the [-1, 1] range,
	// which are stored in the 'points' array. 'h' is the normalized position where a ray hits the fiber.

	// The integration sample points, 'points[i]', represent the normalized hit location 'h' on the fiber.
	// Cache the angle of incidence (gamma_i) for each sample point h.
	std::vector<Float> gammaIs(numPoints);
	for (int i = 0; i < numPoints; ++i) {
		// Convert position h to angle gamma_i via the relation: h = sin(gamma_i).
		gammaIs[i] = std::asin(points[i]);
	}

	// Precompute the Gaussian detector and sample it into three 1D tables.
	// This is the only part of the precomputation that is actually approximate.
	// 2048 samples are enough to support the lowest roughness that the BCSDF
	// can reliably simulate
	const int numGaussianSamples = 2048;
	std::vector<std::vector<Float>> Ds(3, std::vector<Float>(numGaussianSamples));
	for (int p = 0; p < 3; ++p) {
		for (int i = 0; i < numGaussianSamples; ++i) {
			// FIXME: use different beta for different p
			// Generate a sample angle 'phi' in the range [0, 2*PI].
			// The term 'i / (numGaussianSamples - 1.0f)' creates a normalized value from 0.0 to 1.0.
			// Multiplying by 2*PI scales it to the full circle.
			Ds[p][i] = D(_betaR, i / (numGaussianSamples - 1.0f) * 2 * M_PI);
		}
	}

	// Simple wrapped linear interpolation of the precomputed table
	auto approxD = [&](int p, Float phi) {
		// 1. Map the input angle 'phi' (in radians) to a floating-point array index 'u'.
		// The scaling factor (numGaussianSamples - 1) / (2 * PI) converts radians to index units.
		// std::abs handles negative phi, leveraging the symmetry of the D function.
		Float u = std::abs(phi * (0.5f * INV_PI * (numGaussianSamples - 1)));

		// 2. Find the two integer indices that bracket 'u'.
		// x0 is the index of the sample point to the left of our target.
		int x0 = int(u);
		// x1 is the index of the sample point to the right.
		int x1 = x0 + 1;

		// 3. Extract the fractional part of 'u' to use as the interpolation weight.
		// If u was 125.7, it now becomes 0.7.
		u -= x0;

		// 4. Perform linear interpolation between the two sample points.
		// The modulo operator '%' handles wrapping for periodic boundary conditions.
		// correctly interpolating between the last and first elements.
		return Ds[p][x0 % numGaussianSamples] * (1.0f - u) + Ds[p][x1 % numGaussianSamples] * u;
	};

	// Here follows the actual precomputation of the azimuthal scattering functions
	// The scattering functions are parametrized with the azimuthal angle phi,
	// and the cosine of the difference angle, cos(thetaD).
	// This parametrization makes the azimuthal function relatively smooth and allows using
	// really low resolutions for the table (64x64 in this case) without any visual
	// deviation from ground truth, even at the lowest supported roughness setting
	for (int y = 0; y < Resolution; ++y) {
		// Generate a normalized value [0, 1] for the current table row.
		// This value represents cos(difference angle), cos(theta_d), where theta_d = (theta_i - theta_o) / 2.
		Float cosTheta_d = y / (Resolution - 1.0f);

		// Precompute reflection Fresnel factor and reduced absorption coefficient
		// See 5. Azimuthal Scattering
		Float iorPrime = std::sqrt(Eta * Eta - (1.0f - cosTheta_d * cosTheta_d)) / cosTheta_d;

		// The cosine of the transmission angle in the longitudinal plane.
		// FIXME: should it use theta_i and eta?
		Float cosThetaT = std::sqrt(1.0f - (1.0f - cosTheta_d * cosTheta_d) * math::sqr(1.0f / Eta));
		// correspond to mu_a_prime = mu_a / cosTheta after Equation (14)
		// The reduced absorption coefficient
		Vector3 sigma_a_prime = _sigmaA / cosThetaT;

		// Precompute gammaT, f_t and internal absorption across all integration points
		std::vector<Float> fresnelTerms(numPoints);
		std::vector<Float> gammaTs(numPoints);
		std::vector<Vector3> absorptions(numPoints);
		// Similar to precomputing gammaIs, i represents 'h'
		for (int i = 0; i < numPoints; ++i) {
			gammaTs[i] = std::asin(math::clamp(points[i] / iorPrime, -1.0f, 1.0f));
			// Equation (14)
			fresnelTerms[i] = fresnel_dielectric(cosTheta_d * std::cos(gammaIs[i]), 1.0f / Eta);
			absorptions[i] = Vector3(std::exp(-2.0f * sigma_a_prime.x * (1 + std::cos(2 * gammaTs[i]))),
									 std::exp(-2.0f * sigma_a_prime.y * (1 + std::cos(2 * gammaTs[i]))),
									 std::exp(-2.0f * sigma_a_prime.z * (1 + std::cos(2 * gammaTs[i]))));
		}

		for (int phiI = 0; phiI < Resolution; ++phiI) {
			Float phi = 2 * M_PI * phiI / (Resolution - 1.0f);

			Float integralR = 0.0f;
			Vector3 integralTT = Vector3(0.f, 0.f, 0.f);
			Vector3 integralTRT = Vector3(0.f, 0.f, 0.f);

			// Here follows the integration across the fiber width, h.
			// Since we were able to precompute most of the factors that
			// are constant w.r.t phi for a given h,
			// we don't have to do much work here.
			for (int i = 0; i < numPoints; ++i) {
				Float fR = fresnelTerms[i];
				Vector3 T = absorptions[i];

				// Euqation (13)
				Float A_R = fR;
				Vector3 A_TT = (1.0f - fR) * (1.0f - fR) * T;
				Vector3 A_TRT = Vector3(A_TT.x * fR * T.x, A_TT.y * fR * T.y, A_TT.z * fR * T.z);

				integralR += weights[i] * A_R * approxD(0, phi - Phi(gammaIs[i], gammaTs[i], 0));
				integralTT += weights[i] * A_TT * approxD(1, phi - Phi(gammaIs[i], gammaTs[i], 1));
				integralTRT += weights[i] * A_TRT * approxD(2, phi - Phi(gammaIs[i], gammaTs[i], 2));
			}

			valuesR[phiI + y * Resolution] = Vector3(0.5f * integralR, 0.5f * integralR, 0.5f * integralR);
			valuesTT[phiI + y * Resolution] = 0.5f * integralTT;
			valuesTRT[phiI + y * Resolution] = 0.5f * integralTRT;
		}
	}

	// Hand the values off to the helper class to construct sampling CDFs and so forth
	_NR = new PrecomputedAzimuthalLobe(std::move(valuesR));
	_NTT = new PrecomputedAzimuthalLobe(std::move(valuesTT));
	_NTRT = new PrecomputedAzimuthalLobe(std::move(valuesTRT));
}

Vector3f dEonPhaseFunction::anglesToVector(Float theta, Float phi) const {
	Float cosTheta = std::cos(theta);
	return Vector3f(
		std::sin(theta),
		cosTheta * std::cos(phi),
		cosTheta * std::sin(phi));
}

PFPrecalcSample dEonPhaseFunction::integrateBSDF(const Vector3f &wi) const {
	PFPrecalcSample accumulated_value;

	auto pluginManager = PluginManager::getInstance();
	Properties samplerProperties("ldsampler");
	samplerProperties.setInteger("sampleCount", this->m_intergration_sample_cnt);
	ref<ConfigurableObject> samplerObject = pluginManager->createObject(MTS_CLASS(Sampler), samplerProperties);
	ref<Sampler> sampler = dynamic_cast<Sampler *>(samplerObject.get());
	sampler->generate(Point2i(-1));
	if (sampler == nullptr) {
		Log(EError, "Sampler can not be created.");
	}
	for (int i = 1; i <= m_intergration_sample_cnt; i++) {
		Vector3 wo;
		Float pdf;
		auto sampledValue = sampleWithLocalWiForIntegration(wi, wo, pdf, sampler);
		if (!std::isfinite(sampledValue.m_normalizationFactor)) {
			Log(EWarn, "dEonPhaseFunction::integrateBSDF(): skipping non-finite integration sample.");
			sampler->advance();
			continue;
		}
		sampler->advance();
		accumulated_value += sampledValue;
	}
	return accumulated_value / m_intergration_sample_cnt;
}

void dEonPhaseFunction::generateIntegratedHairBSDFTable() {
	Log(EInfo, "Starting to generate integrated BSDF table...");

	if (m_theta_i_steps <= 1) {
		Log(EError, "m_theta_i_steps must be greater than 1, but is %d.", m_theta_i_steps);
		return;
	}

	// Prepare the storage vector.
	m_normalizationFactors.clear();
	m_normalizationFactors.resize(m_theta_i_steps);
	m_totalScatteredEnergy->clear();
	m_totalScatteredEnergy->resize(m_theta_i_steps);
	m_totalTransmissiveEnergy->clear();
	m_totalTransmissiveEnergy->resize(m_theta_i_steps);

	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;

	// Generate table data by integrating for various incident theta angles.
	for (int i = 0; i < m_theta_i_steps; ++i) {
		Float theta_i = min_theta_i + (Float)i / (Float)(m_theta_i_steps - 1) * (max_theta_i - min_theta_i);

		// Assume azimuthal symmetry, so phi_i can be fixed at 0.
		Vector3f wi = anglesToVector(theta_i, 0.0f);

		auto integral_value = integrateBSDF(wi);

		m_normalizationFactors[i] = integral_value.m_normalizationFactor;
		(*m_totalScatteredEnergy)[i] = integral_value.m_totalScatteredEnergy;
		(*m_totalTransmissiveEnergy)[i] = integral_value.m_transmissiveEnergy;

		// Progress reporting.
		if ((i + 1) % (m_theta_i_steps / 10) == 0 || i == m_theta_i_steps - 1) {
			Log(EInfo, "Progress: %d%%", (int)((Float)(i + 1) / m_theta_i_steps * 100.0f));
		}
	}

	Log(EInfo, "Integrated BSDF table generation complete. Table size: %zu.", m_normalizationFactors.size());
}

void dEonPhaseFunction::deriveNormalizationFactor(std::shared_ptr<std::vector<Vector3>> totalScatteredEnergy) {
	m_totalScatteredEnergy = totalScatteredEnergy;
	m_normalizationFactors.clear();
	m_normalizationFactors.resize(m_theta_i_steps);
	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;
	for (int i = 0; i < m_theta_i_steps; ++i) {
		m_normalizationFactors[i] = CwiseMul(m_totalScatteredEnergy->at(i), rgbMask).sum() / maskCount;
	}
}

Float dEonPhaseFunction::lookupNormalizationFactor(Float thetaI) const {
	// Guard against missing precomputed data.
	if (m_normalizationFactors.empty()) {
		Log(EError, "lookupNormalizationFactor called, but the integration table is empty. Returning 0.");
		return 0.0f;
	}

	// Table domain in theta_i.
	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;
	const int table_size = m_normalizationFactors.size();

	// Map theta_i to normalized coordinate u in [0, 1].
	Float u = (thetaI - min_theta_i) / (max_theta_i - min_theta_i);
	u = std::max(0.0f, std::min(1.0f, u)); // Clamp to handle out-of-range inputs

	// Compute interpolation position in table space.
	const Float pos = u * (table_size - 1);

	// Neighboring indices for linear interpolation.
	const int idx0 = static_cast<int>(std::floor(pos));
	const int idx1 = std::min(idx0 + 1, table_size - 1); // Clamp to avoid going out of bounds

	// Interpolation weight (fractional part of pos).
	const Float t = pos - idx0;

	// Values at neighboring indices.
	const Float val0 = m_normalizationFactors[idx0];
	const Float val1 = m_normalizationFactors[idx1];

	// Linear interpolation: (1 - t) * val0 + t * val1
	return val0 + t * (val1 - val0);
}

Vector3 dEonPhaseFunction::lookupTotalScatteredEnergy(Float thetaI) const {
	// Guard against missing precomputed data.
	if (m_totalScatteredEnergy->empty()) {
		Log(EWarn, "lookupNormalizationFactor called, but the integration table is empty. Returning 0.");
		return Vector3(0.0f);
	}

	// Table domain in theta_i.
	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;
	const int table_size = m_totalScatteredEnergy->size();

	// Map theta_i to normalized coordinate u in [0, 1].
	Float u = (thetaI - min_theta_i) / (max_theta_i - min_theta_i);
	u = std::max(0.0f, std::min(1.0f, u)); // Clamp to handle out-of-range inputs

	// Compute interpolation position in table space.
	const Float pos = u * (table_size - 1);

	// Neighboring indices for linear interpolation.
	const int idx0 = static_cast<int>(std::floor(pos));
	const int idx1 = std::min(idx0 + 1, table_size - 1); // Clamp to avoid going out of bounds

	// Interpolation weight (fractional part of pos).
	const Float t = pos - idx0;

	// Values at neighboring indices.
	const Vector3 val0 = (*m_totalScatteredEnergy)[idx0];
	const Vector3 val1 = (*m_totalScatteredEnergy)[idx1];

	// Linear interpolation: (1 - t) * val0 + t * val1
	Vector3 scattered = val0 + t * (val1 - val0);
	scattered = scattered.saturate();
	Vector3 maskedScattered(0.f);
	if (maskAlbedo) {
		maskedScattered = CwiseMul(scattered, rgbMask);
	} else {
		maskedScattered = scattered;
	}

	return maskedScattered;
}

Vector3 dEonPhaseFunction::lookupTotalTransmissiveEnergy(Float thetaI) const {
	// Guard against missing precomputed data.
	if (m_totalTransmissiveEnergy->empty()) {
		Log(EWarn, "lookupNormalizationFactor called, but the integration table is empty. Returning 0.");
		return Vector3(0.0f);
	}

	// Table domain in theta_i.
	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;
	const int table_size = m_totalTransmissiveEnergy->size();

	// Map theta_i to normalized coordinate u in [0, 1].
	Float u = (thetaI - min_theta_i) / (max_theta_i - min_theta_i);
	u = std::max(0.0f, std::min(1.0f, u)); // Clamp to handle out-of-range inputs

	// Compute interpolation position in table space.
	const Float pos = u * (table_size - 1);

	// Neighboring indices for linear interpolation.
	const int idx0 = static_cast<int>(std::floor(pos));
	const int idx1 = std::min(idx0 + 1, table_size - 1); // Clamp to avoid going out of bounds

	// Interpolation weight (fractional part of pos).
	const Float t = pos - idx0;

	// Values at neighboring indices.
	const Vector3 val0 = (*m_totalTransmissiveEnergy)[idx0];
	const Vector3 val1 = (*m_totalTransmissiveEnergy)[idx1];

	// Linear interpolation: (1 - t) * val0 + t * val1
	Vector3 scattered = val0 + t * (val1 - val0);
	scattered = scattered.saturate();
	Vector3 maskedScattered(0.f);
	if (maskAlbedo) {
		maskedScattered = CwiseMul(scattered, rgbMask);
	} else {
		maskedScattered = scattered;
	}

	return maskedScattered;
}

std::string dEonPhaseFunction::generateCacheFilename() const {
	// We combine critical parameters into a string and hash it to keep the filename short
	// but specific to this exact setup.
	std::stringstream ss;
	ss << std::fixed << std::setprecision(4);
	ss << "v1_"; // Versioning prefix
	ss << "_r" << _roughness;
	ss << "_sa" << _scaleAngleDeg;
	ss << "_sig" << _sigmaA.x << "_" << _sigmaA.y << "_" << _sigmaA.z;
	ss << "_steps" << m_theta_i_steps;

	std::hash<std::string> hasher;
	size_t hash = hasher(ss.str());

	std::stringstream filename;
	filename << "dEon_cache_" << std::hex << hash << ".bin";
	return filename.str();
}

bool dEonPhaseFunction::loadPrecomputedTable(const fs::path &path) {
	std::ifstream file(path.string(), std::ios::binary);

	if (!file.is_open()) {
		return false;
	}

	// 1. Read and Check Magic Number
	int magic;
	file.read(reinterpret_cast<char *>(&magic), sizeof(int));

	// Cache format version check.
	if (magic != 0xDE000002) {
		Log(EWarn, "dEonPhaseFunction: Cache file version mismatch (Expected 0xDE000002). Ignoring %s", path.string().c_str());
		return false;
	}

	// 2. Read Parameters and Validate
	Float cachedRoughness, cachedScaleAngle;
	Vector3 cachedSigmaA;
	int cachedSteps;

	file.read(reinterpret_cast<char *>(&cachedRoughness), sizeof(Float));
	file.read(reinterpret_cast<char *>(&cachedScaleAngle), sizeof(Float));
	file.read(reinterpret_cast<char *>(&cachedSigmaA), sizeof(Vector3));
	file.read(reinterpret_cast<char *>(&cachedSteps), sizeof(int));

	if (!math::isApproxEqual(cachedRoughness, _roughness) ||
		!math::isApproxEqual(cachedScaleAngle, _scaleAngleDeg) ||
		!cachedSigmaA.isApproxEqual(_sigmaA) ||
		cachedSteps != m_theta_i_steps) {
		Log(EWarn, "dEonPhaseFunction: Cache file parameters do not match current settings. Ignoring %s", path.string().c_str());
		return false;
	}

	// 3. Read Data

	// Part A: Scattered Energy
	// Ensure the shared pointer is allocated before resizing.
	if (!m_totalScatteredEnergy)
		m_totalScatteredEnergy = std::make_shared<std::vector<Vector3>>();

	m_totalScatteredEnergy->resize(m_theta_i_steps);
	file.read(reinterpret_cast<char *>(m_totalScatteredEnergy->data()), m_theta_i_steps * sizeof(Vector3));

	// Part B: transmissive energy.
	// Ensure the shared pointer is allocated before resizing.
	if (!m_totalTransmissiveEnergy) {
		m_totalTransmissiveEnergy = std::make_shared<std::vector<Vector3>>();
	}

	m_totalTransmissiveEnergy->resize(m_theta_i_steps);
	file.read(reinterpret_cast<char *>(m_totalTransmissiveEnergy->data()), m_theta_i_steps * sizeof(Vector3));

	// Check strict integrity
	if (!file) {
		Log(EError, "dEonPhaseFunction: Error reading data from cache file (file truncated?).");
		return false;
	}

	Log(EInfo, "dEonPhaseFunction: Loaded precomputed table from %s", path.string().c_str());
	return true;
}

void dEonPhaseFunction::savePrecomputedTable(const fs::path &path) const {
	std::ofstream file(path.string(), std::ios::binary);

	if (!file.is_open()) {
		Log(EWarn, "dEonPhaseFunction: Could not save cache file to %s", path.string().c_str());
		return;
	}

	// 1. Write Magic Number/Version
	// Cache format version tag (updated to avoid old-format collisions).
	const int magic = 0xDE000002;
	file.write(reinterpret_cast<const char *>(&magic), sizeof(int));

	// 2. Write Parameters
	file.write(reinterpret_cast<const char *>(&_roughness), sizeof(Float));
	file.write(reinterpret_cast<const char *>(&_scaleAngleDeg), sizeof(Float));
	file.write(reinterpret_cast<const char *>(&_sigmaA), sizeof(Vector3));
	file.write(reinterpret_cast<const char *>(&m_theta_i_steps), sizeof(int));

	// 3. Write Data
	// Validate both vectors before writing.
	size_t size = m_totalScatteredEnergy ? m_totalScatteredEnergy->size() : 0;
	size_t sizeTrans = m_totalTransmissiveEnergy ? m_totalTransmissiveEnergy->size() : 0;

	if (size != (size_t)m_theta_i_steps || sizeTrans != (size_t)m_theta_i_steps) {
		Log(EError, "dEonPhaseFunction: Internal state mismatch during save (vectors empty or wrong size).");
		return;
	}

	// Write Scattered Energy
	file.write(reinterpret_cast<const char *>(m_totalScatteredEnergy->data()), size * sizeof(Vector3));

	// Write transmissive energy.
	file.write(reinterpret_cast<const char *>(m_totalTransmissiveEnergy->data()), sizeTrans * sizeof(Vector3));

	Log(EInfo, "dEonPhaseFunction: Saved precomputed table to %s", path.string().c_str());
	file.close();
}

Float dEonPhaseFunction::integrateBSDFTest(const Vector3f &wi) const {
	auto pluginManager = PluginManager::getInstance();
	Properties samplerProperties("ldsampler");
	samplerProperties.setInteger("sampleCount", this->m_intergration_sample_cnt);
	ref<ConfigurableObject> samplerObject = pluginManager->createObject(MTS_CLASS(Sampler), samplerProperties);
	ref<Sampler> sampler = dynamic_cast<Sampler *>(samplerObject.get());
	if (sampler == nullptr) {
		Log(EError, "Sampler can not be created.");
	}
	sampler->generate(Point2i(-1));

	Float accumulated_value = 0.f;
	Float sinThetaI = wi.x;
	Float thetaI = math::safe_asin(sinThetaI);
	for (int i = 1; i <= m_intergration_sample_cnt; i++) {
		Vector3 wo;
		// Float pdf;
		wo = warp::squareToUniformSphere(sampler->next2D());
		auto sampledValue = evalWithLocalWiWo(wi, wo, 0.f);
		sampledValue = sampledValue / warp::squareToUniformSpherePdf();
		accumulated_value += sampledValue;
		sampler->advance();
	}
	Float C = lookupNormalizationFactor(thetaI);
	if (std::abs(C) < 1e-6 || !std::isfinite(C)) {
		Log(EWarn, "dEonPhaseFunction::integrateBSDFUniform(): invalid normalization factor.");
	}
	accumulated_value = accumulated_value / m_intergration_sample_cnt;
	if (!std::isfinite(accumulated_value)) {
		Log(EWarn, "dEonPhaseFunction::integrateBSDFUniform(): non-finite accumulated value; returning zero.");
		return 0.f;
	}
	return accumulated_value;
}

Object *__dEonPhaseFunction_unSer(Stream *stream, InstanceManager *manager) {
	return new dEonPhaseFunction(stream, manager);
}

DEON_API Class *dEonPhaseFunction::m_theClass = new Class("dEonPhaseFunction", false, "PhaseFunction", NULL, (void *)&__dEonPhaseFunction_unSer);

const Class *dEonPhaseFunction::getClass() const {
	return m_theClass;
}

MTS_EXPORT_PLUGIN(dEonPhaseFunction, "dEon phase function");
MTS_NAMESPACE_END
