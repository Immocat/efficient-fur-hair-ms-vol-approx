#include "dEon.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mitsuba/core/quad.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/shape.h>
MTS_NAMESPACE_BEGIN

// Unimplemented BSDF constructor and interface methods.
// These would typically handle parameter loading from the scene file
// and setting up the derived parameters like _betaR, _vR, etc.

dEonHairBSDF::dEonHairBSDF(const Properties &props)
	: BSDF(props) {
	_roughness = math::clamp(props.getFloat("roughness", 0.5), 0.f, 1.f);
	_melaninConcentration = props.getFloat("melaninConcentration", 0.5);
	_melaninRatio = props.getFloat("melaninRatio", 0.5);
	_scaleAngleDeg = props.getFloat("scaleAngleDeg", 0.f);
	generateIntegrationTable = props.getBoolean("generateIntegrationTable", false);
	useUniformSphereSampling = props.getBoolean("useUniformSphereSampling", false);
	m_maskMode = props.getBoolean("maskMode", false);

	m_overrideSigmaA = props.getBoolean("overrideSigmaA", false);
	if (m_overrideSigmaA) {
		_sigmaA = props.getVector("sigmaA"); // no default value
	}

	calculateDrivedProperty();

	precomputeAzimuthalDistributions();
	precomputeDS();

	if (generateIntegrationTable) {
		std::ostringstream oss;
		oss << "dEon_Rgh-" << _roughness
			<< "_MC-" << std::fixed << std::setprecision(2) << _melaninConcentration
			<< "_MR-" << std::fixed << std::setprecision(2) << _melaninRatio << ".txt";

		std::string fileName = oss.str();

		generateIntegratedHairBSDFTable(
			this,
			fileName,
			90,
			90,
			180);
	}

	configure();
	// generateHairBSDFTable(this, "eval_table", 90, 180);
	// generateHairPDFTable(this, "pdf_table", 90, 180);
	// drive_mitsuba_sample_test(this, 1000, "sample_test");
}

dEonHairBSDF::dEonHairBSDF(Stream *stream, InstanceManager *manager)
	: BSDF(stream, manager) {
	// set default value first
	_roughness = 0.5;
	_melaninConcentration = 0.5;
	_melaninRatio = 0.5;
	calculateDrivedProperty();

	if (!stream->canRead()) {
		return;
	}
	_roughness = math::clamp(stream->readFloat(), 0.f, 1.f);
	_melaninConcentration = stream->readFloat();
	_melaninRatio = stream->readFloat();
	calculateDrivedProperty();
	configure();
}

MatrixRGB dEonHairBSDF::PrecomputeHairMatrixRGB(int bands, int nSamples) const {
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

	// --- Progress Reporting Setup ---
	// Calculate the interval to report (every 10%)
	// Ensure we don't divide by zero if nSamples is small
	int reportInterval = std::max(1, nSamples / 10);

	SLog(EInfo, "Starting Hair Matrix Precomputation (%d samples, %d bands)...", nSamples, bands);

	Intersection its;

	// 3. Monte Carlo Integration Loop
	for (int s = 0; s < nSamples; ++s) {
		// A. Sample Directions
		Vector wi = warp::squareToUniformSphere(sampler->next2D());

		BSDFSamplingRecord bRec(its, sampler.get(), ERadiance);
		bRec.wi = wi;

		Float pdf = 0.f;
		Point2 samplePoint = sampler->next2D();

		// B. Evaluate BSDF (RGB Version)
		Spectrum weight_spec = this->sample(bRec, pdf, samplePoint);

		if(weight_spec.isZero() || pdf<=0.f){
			if((s + 1) % reportInterval == 0){
				Float progress = 100.0f * static_cast<Float>(s + 1) / static_cast<Float>(nSamples);
				SLog(EInfo, "Precomputing... %.0f%% complete", progress);
			}
			sampler->advance();
			continue;
		}

		Vector wo = bRec.wo;

		Float w_r, w_g, w_b;

		weight_spec.toLinearRGB(w_r, w_g, w_b);

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

		// --- Progress Reporting ---
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

void dEonHairBSDF::configure() {
	m_components.clear();
	// Register this BSDF as a glossy reflection model
	// This tells the integrator that it can be sampled and evaluated.
	int flags = EDiffuseReflection | EUsesSampler | EFrontSide;
	m_components.push_back(flags);

	// Call the parent class's configure last
	BSDF::configure();
}

Spectrum dEonHairBSDF::sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
	Float pdf;
	Spectrum value = this->sample(bRec, pdf, sample);
	return value;
}

Spectrum dEonHairBSDF::sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
	if (useUniformSphereSampling) {
		bRec.wo = warp::squareToUniformSphere(sample);
		pdf = warp::squareToUniformSpherePdf();
		bRec.eta = Eta;
		bRec.sampledType = EGlossyReflection;
		return eval(bRec, ESolidAngle) / pdf;
	}

	Point2 xiN = sample;
	Point2 xiM = bRec.sampler->next2D();

	const Vector3f &wi = bRec.wi;

	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = math::safe_asin(sinThetaI);

	Float phiI = std::atan2(wi.z, wi.y);

	//  Lobe shift due to hair scale tilt, following the values in
	//  "Importance Sampling for Physically-Based Hair Fiber Models"
	//  rather than the earlier paper by Marschner et al. I believe
	//  these are slightly more accurate.
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// The following lines are just lobe selection
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
	Float sinPhiO = std::sin(phi + phiI);
	Float cosPhiO = std::cos(phi + phiI);

	Vector3 wo = Vector3(
		sinThetaO,
		cosPhiO * cosThetaO,
		sinPhiO * cosThetaO);

	bRec.wo = wo;
	bRec.eta = Eta;
	bRec.sampledType = EGlossyReflection;

	pdf = this->pdf(bRec, EMeasure::ESolidAngle);
	if (!std::isfinite(pdf) || std::isnan(pdf)) {
		Log(EWarn, "dEonHairBSDF::sample(): non-finite PDF; returning zero contribution.");
		pdf = 0.f;
		return Spectrum(0.f);
	}
	if (pdf > 0.f) {
		Spectrum result = this->eval(bRec, EMeasure::ESolidAngle) / pdf;
		Vector t;
		result.toLinearRGB(t.x, t.y, t.z);
		if (!std::isfinite(t.x) || std::isnan(t.x)) {
			Log(EWarn, "dEonHairBSDF::sample(): non-finite sample value; returning zero contribution.");
			return Spectrum(0.f);
		}
		return result;
	} else {
		return Spectrum(0.f);
	}
}

Spectrum dEonHairBSDF::sampleWithWorldWiAndOrientation(const Vector3 &worldWi, const Vector3 &orientation, Float &pdf, Vector3f &worldWo, Sampler *sampler) const {
	if (orientation.isZero()) {
		Log(EError, "Orientation is zero when sampling Hair BSDF!");
	}
	Frame frame = Frame::fromS(orientation);
	Vector wi = frame.toLocal(worldWi);

	Vector wo_local;
	Float _eta;
	BSDF::EBSDFType _sampledType;
	Spectrum sampleValue = sampleWithLocalWi(wi, wo_local, pdf, sampler->next2D(), _eta, _sampledType, sampler);
	worldWo = frame.toWorld(wo_local);
	return sampleValue;
}

Spectrum dEonHairBSDF::sampleWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, const Point2 &sample, Float &eta, BSDF::EBSDFType &sampledType, Sampler *sampler) const {
	Point2 xiN = sample;
	Point2 xiM = sampler->next2D();

	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = math::safe_asin(sinThetaI);

	Float phiI = std::atan2(wi.z, wi.y);

	//  Lobe shift due to hair scale tilt, following the values in
	//  "Importance Sampling for Physically-Based Hair Fiber Models"
	//  rather than the earlier paper by Marschner et al. I believe
	//  these are slightly more accurate.
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// The following lines are just lobe selection
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
	Float sinPhiO = std::sin(phi + phiI);
	Float cosPhiO = std::cos(phi + phiI);
	wo = Vector3(
		sinThetaO,
		cosPhiO * cosThetaO,
		sinPhiO * cosThetaO);

	eta = Eta;
	sampledType = EGlossyReflection;

	pdf = this->pdfWithLocalWiWo(wi, wo);
	if (!std::isfinite(pdf) || std::isnan(pdf)) {
		Log(EWarn, "dEonHairBSDF::sampleWithLocalWi(): non-finite PDF; returning zero contribution.");
		pdf = 0.f;
		return Spectrum(0.f);
	}
	if (pdf > 0.f) {
		Spectrum result = this->evalWithLocalWiWo(wi, wo) / pdf;
		Vector t;
		result.toLinearRGB(t.x, t.y, t.z);
		if (!std::isfinite(t.x) || std::isnan(t.x)) {
			Log(EWarn, "dEonHairBSDF::sampleWithLocalWi(): non-finite sample value; returning zero contribution.");
			return Spectrum(0.f);
		}
		return result;
	} else {
		return Spectrum(0.f);
	}
}

Float dEonHairBSDF::pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
	if (useUniformSphereSampling) {
		return warp::squareToUniformSpherePdf();
	}
	return pdfWithLocalWiWo(bRec.wi, bRec.wo);
}

Float dEonHairBSDF::pdfWithLocalWiWo(const Vector3 &wi, const Vector3 &wo) const {
	if (useUniformSphereSampling) {
		return warp::squareToUniformSpherePdf();
	}

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

	//  Lobe shift due to hair scale tilt, following the values in
	//  "Importance Sampling for Physically-Based Hair Fiber Models"
	//  rather than the earlier paper by Marschner et al. I believe
	//  these are slightly more accurate.
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

	return (1.0f / weightSum) *
		   (pdfR * _NR->pdf(phi, cosThetaD) + pdfTT * _NTT->pdf(phi, cosThetaD) + pdfTRT * _NTRT->pdf(phi, cosThetaD));
}

Spectrum dEonHairBSDF::eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
	// Spectrum r;
	// auto v = (bRec.its.shFrame.s + Vector3(1.f)) / 2;
	// r.fromLinearRGB(v.x, v.y, v.z);
	//////r.fromLinearRGB(1.f, 1.f, 1.f);
	// return r;
	if (m_maskMode) {
		return Spectrum(1.f);
	}

	const Vector3f &wi = bRec.wi;
	const Vector3f &wo = bRec.wo;

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

	//  Lobe shift due to hair scale tilt, following the values in
	//  "Importance Sampling for Physically-Based Hair Fiber Models"
	//  rather than the earlier paper by Marschner et al. I believe
	//  these are slightly more accurate.
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Evaluate longitudinal scattering functions
	Float MR = M(_vR, std::sin(thetaIR), sinThetaO, std::cos(thetaIR), cosThetaO);
	Float MTT = M(_vTT, std::sin(thetaITT), sinThetaO, std::cos(thetaITT), cosThetaO);
	Float MTRT = M(_vTRT, std::sin(thetaITRT), sinThetaO, std::cos(thetaITRT), cosThetaO);

	// Vector3 NR_eval = _NR->eval(phi, cosThetaD);
	// Vector3 NTT_eval = _NTT->eval(phi, cosThetaD);
	// Vector3 NTRT_eval = _NTRT->eval(phi, cosThetaD);

	Vector3 R = MR * _NR->eval(phi, cosThetaD);
	Vector3 TT = MTT * _NTT->eval(phi, cosThetaD);
	Vector3 TRT = MTRT * _NTRT->eval(phi, cosThetaD);
	Vector3 combinedLobes = R + TT + TRT;
	// if (cosThetaD < 0.3) {
	//	return Spectrum(0.f);
	// }

	Vector3 S = (combinedLobes);
	Spectrum ret;
	ret.fromLinearRGB(S.x, S.y, S.z);

	return ret;
}

Spectrum dEonHairBSDF::evalWithLocalWiWo(const Vector3 &wi, const Vector3 &wo) const {
	// Spectrum r;
	// auto v = (bRec.its.shFrame.s + Vector3(1.f)) / 2;
	// r.fromLinearRGB(v.x, v.y, v.z);
	//////r.fromLinearRGB(1.f, 1.f, 1.f);
	// return r;
	if (m_maskMode) {
		return Spectrum(1.f);
	}

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

	//  Lobe shift due to hair scale tilt, following the values in
	//  "Importance Sampling for Physically-Based Hair Fiber Models"
	//  rather than the earlier paper by Marschner et al. I believe
	//  these are slightly more accurate.
	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	// Evaluate longitudinal scattering functions
	Float MR = M(_vR, std::sin(thetaIR), sinThetaO, std::cos(thetaIR), cosThetaO);
	Float MTT = M(_vTT, std::sin(thetaITT), sinThetaO, std::cos(thetaITT), cosThetaO);
	Float MTRT = M(_vTRT, std::sin(thetaITRT), sinThetaO, std::cos(thetaITRT), cosThetaO);

	// Vector3 NR_eval = _NR->eval(phi, cosThetaD);
	// Vector3 NTT_eval = _NTT->eval(phi, cosThetaD);
	// Vector3 NTRT_eval = _NTRT->eval(phi, cosThetaD);

	Vector3 R = MR * _NR->eval(phi, cosThetaD);
	Vector3 TT = MTT * _NTT->eval(phi, cosThetaD);
	Vector3 TRT = MTRT * _NTRT->eval(phi, cosThetaD);
	Vector3 combinedLobes = R + TT + TRT;
	// if (cosThetaD < 0.3) {
	//	return Spectrum(0.f);
	// }

	Vector3 S = (combinedLobes);
	Spectrum ret;
	ret.fromLinearRGB(S.x, S.y, S.z);

	return ret;
}

// #pragma optimize("", off)
Spectrum dEonHairBSDF::forwardScattering(const Scene &scene, const Float &hairThickness, const BSDFSamplingRecord &bRec, const Vector3 &wd, Float dist, const Point &safeOrigin, Vector3 &sigma_f, Vector3 &trans, bool &direct) const {
	const Vector3f &wi = bRec.wi;

	Float sinThetaD = wd.x;
	Float cosThetaD = math::trigInverse(sinThetaD);
	Float SinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(SinThetaI);
	Float thetaD = std::asin(std::clamp(sinThetaD, -1.f, 1.f));
	Float thetaI = std::asin(std::clamp(SinThetaI, -1.f, 1.f));
	Float theta_d = (thetaD - thetaI) / 2;

	Vector3 sigma_f_sqr(0.f);
	trans = forwardTrans(scene, bRec, wi, wd, hairThickness, dist, safeOrigin, sigma_f_sqr, direct);
	sigma_f = sigma_f_sqr.sqrt();

	Vector3 s_f(0.f);
	// Float theta_h = (thetaD + thetaI) / 2.f;
	//  original Paper's section3/4 implementation
	Float phi_d = std::atan2(wd.z, wd.y);
	Float phi_i = std::atan2(wi.z, wi.y);
	Float phi = phi_d - phi_i;
	if (phi < 0.f) {
		phi += 2.0f * M_PI;
	}

	// The S term in the section 3.1.2, but it is not used in the final implementation (section 4.2)
	/*if (trans == Vector3(1.f)) {
		s_f = Vector3(1.f);
	} else {
		if ((phi > 0.5f * M_PI) && (phi < 1.5f * M_PI)) {
			Float s_f_r = 1 / M_PI / cos(theta_d) * D(sigma_f[0], thetaD + thetaI);
			Float s_f_g = 1 / M_PI / cos(theta_d) * D(sigma_f[1], thetaD + thetaI);
			Float s_f_b = 1 / M_PI / cos(theta_d) * D(sigma_f[2], thetaD + thetaI);
			s_f = Vector3(s_f_r, s_f_g, s_f_b);
		} else {
			s_f = Vector3(0.f);
		}
	}

	trans = CwiseMul(trans, s_f);*/

	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

	Intersection it;
	Spectrum result;
	Vector3 S(0.f);
	if (direct) {
		Float M_R = M(_betaR * _betaR, sin(thetaIR), sinThetaD, cos(thetaIR), cosThetaD);
		Vector3 N_R = _NR->eval(phi, cos(theta_d));
		Float M_TT = M(_betaTT * _betaTT, sin(thetaITT), sinThetaD, cos(thetaITT), cosThetaD);
		Vector3 N_TT = _NTT->eval(phi, cos(theta_d));
		Float M_TRT = M(_betaTRT * _betaTRT, sin(thetaITRT), sinThetaD, cos(thetaITRT), cosThetaD);
		Vector3 N_TRT = _NTRT->eval(phi, cos(theta_d));

		S = (M_R * N_R + M_TT * N_TT + M_TRT * N_TRT);
	} else {
		Float M_R = M(_betaR * _betaR + sigma_f_sqr[0], sin(thetaIR), sinThetaD, cos(thetaIR), cosThetaD);
		Vector3 N_R = _NR->eval(phi, cos(theta_d));
		Float M_TT = M(_betaTT * _betaTT + sigma_f_sqr[0], sin(thetaITT), sinThetaD, cos(thetaITT), cosThetaD);
		Vector3 N_TT = _NTT->eval(phi, cos(theta_d));
		Float M_TRT = M(_betaTRT * _betaTRT + sigma_f_sqr[0], sin(thetaITRT), sinThetaD, cos(thetaITRT), cosThetaD);
		Vector3 N_TRT = _NTRT->eval(phi, cos(theta_d));

		// trans * dense_factor * widen lobe
		S = (M_R * N_R + M_TT * N_TT + M_TRT * N_TRT);
		S = Vector3(
			df * S.x * trans.x,
			df * S.y * trans.y,
			df * S.z * trans.z);
	}

	result.fromLinearRGB(S.x, S.y, S.z);
	return result;
}

Vector3 dEonHairBSDF::forwardTrans(const Scene &scene, const BSDFSamplingRecord &bRec,
								   const Vector3 &wi, const Vector3 &wd,
								   const Float hairThickness, Float dist,
								   const Point &safeOrigin,
								   Vector3 &var, bool &direct) const {
	Point p = bRec.its.p;
	Vector3 trans(1.f, 1.f, 1.f);

	Vector3 world_wd = bRec.its.toWorld(wd);
	Point nextP = p;

	Ray shadow_ray(safeOrigin, world_wd, 0);
	shadow_ray.mint = Epsilon;
	shadow_ray.maxt = dist * (1.0f - ShadowEpsilon);

	var = Vector3(0, 0, 0);
	direct = true;
	int cur = 0;

	Intersection shadowRayIts;

	// Walk through hair intersections along the shadow ray.
	while (scene.rayIntersect(shadow_ray, shadowRayIts)) {
		if (shadowRayIts.isEmitter()) {
			break;
		}

		direct = false;

		bool isHair = false;
		if (shadowRayIts.shape && shadowRayIts.shape->isCurve()) {
			isHair = true;
		}

		if (!isHair) {
			direct = false;
			return Vector3(0.0f);
		}

		Vector3 wi = -shadow_ray.d.normalized();
		Vector3 wi_local = shadowRayIts.toLocal(wi);

		Float sinThetaI = wi_local.x;
		Float theta_i = std::asin(std::clamp(sinThetaI, -1.0f, 1.0f));
		Float theta_i_bin_Float = (theta_i + M_PI / 2.0f) / M_PI * NUM_BINS - 0.5f;
		int k = std::clamp((int)theta_i_bin_Float, 0, NUM_BINS - 1);

		Vector3 trans_cur = a_f[k];
		trans = CwiseMul(trans, trans_cur);
		var += Vector3(beta_f_sqr[k]);

		if (trans.max() < 1e-3f)
			return Vector3(0.0f);

		Vector3 local_wd = shadowRayIts.toLocal(world_wd);
		if (dot(world_wd, shadowRayIts.shFrame.n) < 0) {
			nextP = shadowRayIts.p - shadowRayIts.shFrame.n.normalized() * (hairThickness);
		} else {
			nextP = shadowRayIts.p + shadowRayIts.shFrame.n.normalized() * (hairThickness * 0.1f);
		}

		shadow_ray = Ray(nextP, world_wd, shadow_ray.time);
		shadow_ray.mint = Epsilon;

		Float distToLight = dist - (nextP - p).length();
		if (distToLight <= 0)
			break;
		shadow_ray.maxt = distToLight * (1.0f - ShadowEpsilon);

		cur++;
		if (cur > 50)
			break;
	}

	if (direct)
		return Vector3(1, 1, 1);
	else {
		return trans * df;
	}
}

Spectrum dEonHairBSDF::backwardScattering(const BSDFSamplingRecord &bRec, const Vector3 &sigmaf, const bool direct, const Vector3 &trans, Vector3 &sigma_b) const {
	Point3 p = bRec.its.p;
	const Vector3 &wi = bRec.wi;
	const Vector3 &wo = bRec.wo;
	Float sinThetaO = std::clamp(wo.x, -1.0f, 1.0f);
	Float sinThetaI = std::clamp(wi.x, -1.0f, 1.0f);
	Float cosThetaO = math::trigInverse(sinThetaO);
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = std::asin(std::clamp(sinThetaI, -1.0f, 1.0f));
	Float thetaO = std::asin(std::clamp(sinThetaO, -1.0f, 1.0f));
	Float theta_h = (thetaO + thetaI) * 0.5f;

	// 3.2.1 Average Backscattering Attenuation
	Float theta_d = (thetaO - thetaI) / 2;
	Float theta_d_bin_Float = (theta_d + M_PI / 2.0f) / M_PI * NUM_BINS - 0.5f;
	int k = (int)floor(theta_d_bin_Float);

	Vector3 avg_a_b = a_b[k];
	Vector3 avg_a_f = a_f[k];
	// Equation (11)
	Vector3 avgA1 = CwiseDiv(CwiseMul(avg_a_b, avg_a_f.pow(2)), (Vector3(1.f) - avg_a_f.pow(2)));
	// Equation (13)
	Vector3 avgA3 = CwiseDiv(CwiseMul(avg_a_b.pow(3), avg_a_f.pow(2)), (Vector3(1.f) - avg_a_f.pow(2)).pow(3));
	// Equation (14)
	Vector3 avgA_B = avgA1 + avgA3;

	// 3.2.2 Average Backscattering Spread
	// avgDelta_b and avg_sigma_b is precomputed
	Float phi_o = std::atan2(wo.z, wo.y);
	Float phi_i = std::atan2(wi.z, wi.y);
	Float phi = phi_o - phi_i;
	if (phi < 0.f)
		phi += 2.0f * M_PI;

	Spectrum ret;
	Vector3 localScatter(0.f);
	Vector3 S_B_Hat(0.f);

	sigma_b = (sigma_b_sqr[k] + sigmaf.sqr()).sqrt();

	// The paper's formula g(...) uses variance, but our implementation takes standard deviation
	if ((phi < 0.5f * M_PI) || (phi > 1.5f * M_PI)) {
		Float S_B_Hat_x = D(sigma_b[0], 2.0f * (theta_h - Delta_b[k][0]));
		Float S_B_Hat_y = D(sigma_b[1], 2.0f * (theta_h - Delta_b[k][1]));
		Float S_B_Hat_z = D(sigma_b[2], 2.0f * (theta_h - Delta_b[k][2]));
		S_B_Hat = Vector3(S_B_Hat_x, S_B_Hat_y, S_B_Hat_z);
	} else {
		S_B_Hat = Vector3(0, 0, 0);
	}
	Float cosThetaD = cos(theta_d);
	if (std::abs(cosThetaD) < Epsilon) {
		return Spectrum(0.f);
	}

	if (direct) {
		localScatter = db * 2.f * CwiseMul(avgA_B, S_B_Hat) / M_PI / cos(theta_d) / cos(theta_d);
	} else {
		localScatter = df * db * 2.f * CwiseMul(avgA_B, S_B_Hat) / M_PI / cos(theta_d) / cos(theta_d);
		localScatter = CwiseMul(trans, localScatter);
	}

	ret.fromLinearRGB(localScatter.x, localScatter.y, localScatter.z);
	return ret;
}
// #pragma optimize("", on)
void dEonHairBSDF::calculateDrivedProperty() {
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

/**
 * Implements Equation (15) from the d'Eon paper.
 * I₀(x) ≈ Σ [ (x²/4)ⁱ / (i!)² ] for i = 0 to 10.
 * This function is a core component of the longitudinal scattering function M_p.
 */
Float dEonHairBSDF::I0(Float x) const {
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

/**
 * Computes log(I₀(x)) with a stable approximation for large x.
 * The standard I0(x) function can overflow for large x (e.g., x > 30).
 * This function switches to an asymptotic expansion of log(I₀(x)) for x > 12
 * to maintain numerical stability, a common trick when dealing with Bessel functions.
 */
Float dEonHairBSDF::logI0(Float x) const {
	if (x > 12.0f) {
		// Asymptotic expansion for log(I₀(x)) for large x:
		// log(I₀(x)) ≈ x - 0.5*log(2πx) + 1/(8x)
		return x + 0.5f * (std::log(1.0f / (2 * M_PI * x)) + 1.0f / (8.0f * x));
	} else {
		return std::log(I0(x));
	}
}

// Standard Gaussian probability density function, g(β; x)
Float dEonHairBSDF::g(Float beta, Float theta) const {
	return std::exp(-theta * theta / (2.0f * beta * beta)) / (math::safe_sqrt(2.0f * M_PI) * beta);
}

/**
 * Implements Equation (11): D_p(φ) = Σ g(β_p; φ - 2πk) for k = -∞ to ∞.
 * This "detector" function models azimuthal blur due to roughness.
 * The infinite sum is approximated by summing terms until they become negligible.
 * It handles the periodic nature of angles by summing over shifted Gaussians.
 */
Float dEonHairBSDF::D(Float beta, Float phi) const {
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

// Implements Equation (8): Φ(p,h) = 2pγ_t - 2γ_i + pπ
Float dEonHairBSDF::Phi(Float gammaI, Float gammaT, int p) const {
	return 2.0f * p * gammaT - 2.0f * gammaI + p * M_PI;
}

// Standard unpolarized Fresnel equation for dielectric materials.
// Used to calculate reflection and transmission at interfaces.
inline Float fresnel_dielectric(Float cos_theta_i, Float eta) {
	Float sin_theta_t_2 = eta * eta * (1.0f - cos_theta_i * cos_theta_i);

	// Total internal reflection
	if (sin_theta_t_2 > 1.0)
		return 1.0;

	Float cos_theta_t = std::sqrt(std::max(1.0f - sin_theta_t_2, 0.0f));

	Float rs = (eta * cos_theta_t - cos_theta_i) / (eta * cos_theta_t + cos_theta_i);
	Float rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

	return 0.5f * (rs * rs + rp * rp);
}

// The following two functions compute the integrands for the azimuthal scattering
// function N_p from Equation (10). As noted by the original author, they are
// primarily for reference or for use in an offline precomputation step. The
// runtime evaluation uses a precomputed lookup table instead.

// Integrand for N₀ (R-lobe, p=0): A(0,h) * D₀(φ - Φ(0,h))
Float dEonHairBSDF::NrIntegrand(Float beta, Float halfWiDotWo, Float phi, Float h) const {
	// γ_i = arcsin(h)
	Float gammaI = std::asin(math::clamp(h, -1.0f, 1.0f));

	// The azimuthal shift is φ - Φ(0,h) = φ - (-2γ_i) = φ + 2γ_i
	Float deltaPhi = phi + 2.0f * gammaI;
	deltaPhi = std::fmod(deltaPhi, 2 * M_PI); // Wrap angle to [0, 2pi)
	if (deltaPhi < 0.0f)
		deltaPhi += 2 * M_PI;

	// A(0,h) from Eq(12) is just the Fresnel term F(η, arccos(ωᵢ⋅ωᵣ)/2)
	Float A0 = fresnel_dielectric(halfWiDotWo, 1.0 / Eta);

	return D(beta, deltaPhi) * A0;
}

// Integrand for Nₚ (internal lobes, p>=1): A(p,h) * Dₚ(φ - Φ(p,h))
Vector3 dEonHairBSDF::NpIntegrand(Float beta, Float cosThetaD, Float phi, int p, Float h) const {
	// η' = sqrt(η² - sin²θ_d) / cosθ_d (effective index of refraction)
	Float iorPrime = std::sqrt(Eta * Eta - (1.0f - cosThetaD * cosThetaD)) / cosThetaD;
	Float cosThetaT = std::sqrt(1.0f - (1.0f - cosThetaD * cosThetaD) * math::sqr(1.0f / Eta));
	// µₐ' = µₐ / cosθ_t (absorption coefficient adjusted for path length)
	Vector3 sigmaAPrime = _sigmaA / cosThetaT;

	// γ_i = arcsin(h)
	Float gammaI = std::asin(math::clamp(h, -1.0f, 1.0f));
	// γ_t = arcsin(h/η')
	Float gammaT = std::asin(math::clamp(h / iorPrime, -1.0f, 1.0f));

	// Internal path length for absorption. The author notes this is a corrected
	// version compared to the one in the original Marschner and d'Eon papers.
	Float l = 2.0f * std::cos(gammaT);

	// f = F(η, arccos(cos(θ_d)cos(arcsin(h)))), internal Fresnel reflectance (Eq 14)
	Float f = fresnel_dielectric(cosThetaD * math::trigInverse(h), 1.0f / Eta);
	// T = exp(-2µₐ' (1 + cos(2γ_t))), but the paper is wrong, should be exp(-µₐ' * l)
	Vector3 T = Vector3(
		std::exp(-sigmaAPrime.x * l),
		std::exp(-sigmaAPrime.y * l),
		std::exp(-sigmaAPrime.z * l));

	// A(p,h) = (1-f)² * f^(p-1) * T^p (Eq 13)
	// Start with (1-f)² * T
	Vector3 Aph = (1.0f - f) * (1.0f - f) * T;
	// Multiply by f and T for the remaining p-1 internal bounces
	for (int i = 1; i < p; ++i) {
		Aph = Vector3(Aph.x * f * T.x, Aph.y * f * T.y, Aph.z * f * T.z);
	}

	// Azimuthal shift φ - Φ(p,h) using the general formula
	Float deltaPhi = phi - Phi(gammaI, gammaT, p);
	deltaPhi = std::fmod(deltaPhi, 2 * M_PI); // Wrap angle
	if (deltaPhi < 0.0f)
		deltaPhi += 2 * M_PI;

	return Aph * D(beta, deltaPhi);
}

/**
 * Implements Equation (7): The longitudinal scattering function M_p.
 * M_p = (csch(1/v) / 2v) * exp( (sin(-θᵢ)sin(θᵣ) / v) ) * I₀( (cos(-θᵢ)cos(θᵣ) / v) )
 * Note: sin(-θᵢ) = -sin(θᵢ) and cos(-θᵢ) = cos(θᵢ).
 */
Float dEonHairBSDF::M(Float v, Float sinThetaI, Float sinThetaO, Float cosThetaI, Float cosThetaO) const {
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

/**
 * Importance samples the M_p distribution to generate a scattered angle θ_o.
 * This is crucial for efficient path tracing of hair. The code uses a more
 * numerically stable inversion of the CDF than the one derived naively from the paper.
 * The technique is from "Numerically stable sampling of the von Mises Fisher distribution...".
 */
Float dEonHairBSDF::sampleM(Float v, Float sinThetaI, Float cosThetaI, Float xi1, Float xi2) const {
	// Sample the polar angle component using the stable inversion formula.
	// This gives us a new `cosTheta` relative to the fiber tangent.
	Float cosTheta = 1.0f + v * std::log(xi1 + (1.0f - xi1) * std::exp(-2.0f / v));
	Float sinTheta = math::trigInverse(cosTheta);

	// Sample the azimuthal angle component uniformly.
	Float cosPhi = std::cos(2 * M_PI * xi2);

	// Combine the sampled components with the incident angle using spherical trigonometry
	// to get the final outgoing longitudinal angle sin(θ_o).
	return -cosTheta * sinThetaI + sinTheta * cosPhi * cosThetaI;
}

void dEonHairBSDF::precomputeAzimuthalDistributions() {
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

inline void regularizeTheta(Float &theta) {
	if (theta > M_PI / 2.0f)
		theta = M_PI - theta;
	else if (theta < -M_PI / 2.0f)
		theta = -M_PI - theta;
}

inline void regularizePhi(Float &phi) {
	while (phi > M_PI)
		phi -= 2.0f * M_PI;
	while (phi < -M_PI)
		phi += 2.0f * M_PI;
}

void dEonHairBSDF::precomputeDS() {

	const int Resolution = PrecomputedAzimuthalLobe::AzimuthalResolution;

	// Iterate over discretized incident longitudinal angles (theta_d).
	// Loop over discretized incident longitudinal angles theta_d.
	for (int k = 0; k < NUM_BINS; k++) {
		// Incident longitudinal angle theta_d, ranging from -PI/2 to PI/2.
		Float theta_d = (k + 0.5f) / NUM_BINS * M_PI - M_PI * 0.5f;

		// The backward lobe ('b') is an average of R and TRT. The forward lobe ('f') is TT.

		// Variance (beta^2) for the backward and forward lobes.
		beta_b_sqr[k] = pow((_betaR + _betaTRT) * 0.5f, 2.0f);
		beta_f_sqr[k] = pow(_betaTT, 2.0f);

		// Longitudinal shifts (alpha) for the primary Marschner lobes.
		Float alpha_R = 2.0f * _scaleAngleRad;
		Float alpha_TT = -_scaleAngleRad;
		Float alpha_TRT = -4.0f * _scaleAngleRad;

		// Mean longitudinal shift for the combined backward and forward lobes.
		alpha_b[k] = (alpha_R + alpha_TRT) / 2.0f;
		alpha_f[k] = alpha_TT;

		// Numerically integrate the average backward scattering attenuation (corresponds to Eq. 12 in the paper) ---
		int NUM_THETA_SAMPLES = NUM_BINS * 2;
		int NUM_PHI_SAMPLES = 12;
		a_b[k] = Vector3(0.0f, 0.0f, 0.0f);
		// Integrate over the outgoing solid angle d_omega (loops i and j).
		for (int i = 0; i < NUM_THETA_SAMPLES; i++) {
			// Outgoing longitudinal angle theta_r.
			Float theta_r = (i + 0.5f) / NUM_THETA_SAMPLES * M_PI - M_PI / 2.0f;
			for (int j = 0; j < NUM_PHI_SAMPLES; j++) {
				// Outgoing azimuthal angle phi_r. The range [-PI/2, PI/2] defines the backward hemisphere.
				Float phi_r = (j + 0.5f) / NUM_PHI_SAMPLES * M_PI - M_PI / 2.0f;
				regularizePhi(phi_r);
				// Differential solid angle d_omega = cos(theta_o) * d_theta_o * d_phi_o. (TODO: is this correct?)
				Float dOmega = cos(theta_r) * (M_PI / NUM_THETA_SAMPLES) * (M_PI / NUM_PHI_SAMPLES);

				// Integrate over the incoming azimuthal angle phi_i (loop l).
				for (int l = 0; l < NUM_PHI_SAMPLES; l++) {
					Float phi_i = (l + 0.5f) / NUM_PHI_SAMPLES * M_PI - M_PI / 2.0f;
					Vector3 wo = Vector3(sin(theta_r), cos(theta_r) * cos(phi_r), cos(theta_r) * sin(phi_r));
					Vector3 wi = Vector3(sin(theta_d), cos(theta_d) * cos(phi_i), cos(theta_d) * sin(phi_i));
					Intersection _;
					BSDFSamplingRecord bRec(_, wi, wo);
					// Evaluate the single-scattering BSDF, f_s().
					Spectrum S = eval(bRec, ESolidAngle);
					Vector3 f(0.f);
					S.toLinearRGB(f.x, f.y, f.z);
					// Accumulate the integrand: f_s * cos(theta_d) * d_phi_i * d_omega.
					a_b[k] += f * cos(theta_d) * (M_PI / NUM_PHI_SAMPLES) * dOmega;
				}
			}
		}
		// Apply the 1/PI normalization factor from the formula.
		a_b[k] /= M_PI;

		// Numerically integratethe average forward scattering attenuation (corresponds to Eq. 6 in the paper) ---
		// This is identical to the a_b calculation, but integrates over the forward hemisphere.
		a_f[k] = Vector3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < NUM_THETA_SAMPLES; i++) {
			Float theta_r = (i + 0.5f) / NUM_THETA_SAMPLES * M_PI - M_PI / 2.0f;
			for (int j = 0; j < NUM_PHI_SAMPLES; j++) {
				// Outgoing azimuthal angle phi_o. The range [PI/2, 3*PI/2] defines the forward hemisphere.
				Float phi_r = (j + 0.5f) / NUM_PHI_SAMPLES * M_PI + M_PI / 2.0f;
				regularizePhi(phi_r);
				Float dOmega = cos(theta_r) * (M_PI / NUM_THETA_SAMPLES) * (M_PI / NUM_PHI_SAMPLES);
				for (int l = 0; l < NUM_PHI_SAMPLES; l++) {
					Float phi_i = (l + 0.5f) / NUM_PHI_SAMPLES * M_PI - M_PI / 2.0f;
					Vector3 wo = Vector3(sin(theta_r), cos(theta_r) * cos(phi_r), cos(theta_r) * sin(phi_r));
					Vector3 wi = Vector3(sin(theta_d), cos(theta_d) * cos(phi_i), cos(theta_d) * sin(phi_i));
					Intersection _;
					BSDFSamplingRecord bRec(_, wi, wo);
					Spectrum S = eval(bRec, ESolidAngle);
					Vector3 f(0.f);
					S.toLinearRGB(f.x, f.y, f.z);
					a_f[k] += f * cos(theta_d) * (M_PI / NUM_PHI_SAMPLES) * dOmega;
				}
			}
		}
		a_f[k] /= M_PI;

		// Calculate A_b. Equation (11) (13) (14)
		Vector3 A_1, A_3;
		for (int i = 0; i < 3; i++) {
			A_1[i] = a_b[k][i] * pow(a_f[k][i], 2.0f) / (1.0f - pow(a_f[k][i], 2.0f));
			A_3[i] = pow(a_b[k][i], 3.0f) * pow(a_f[k][i], 2.0f) / pow(1.0f - pow(a_f[k][i], 2.0f), 3.0f);
		}
		A_b[k] = A_1 + A_3;

		// Calculate sigma_b. Equation (17)
		Vector3 sigma_b;
		for (int i = 0; i < 3; i++) {
			sigma_b[i] = (1.0f + db * a_f[k][i] * a_f[k][i]);
			sigma_b[i] *= (a_b[k][i] * sqrt(2.0f * beta_f_sqr[k] + beta_b_sqr[k]) +
						   a_b[k][i] * a_b[k][i] * a_b[k][i] * sqrt(2.0f * beta_f_sqr[k] + 3.0f * beta_b_sqr[k])) /
						  (a_b[k][i] + a_b[k][i] * a_b[k][i] * a_b[k][i] *
										   (2.0f * sqrt(beta_f_sqr[k]) + 3.0f * sqrt(beta_b_sqr[k])));
			sigma_b_sqr[k][i] = sigma_b[i] * sigma_b[i];
		}

		// Calculate Delta_b.
		Delta_b[k] = Vector3(0.0f, 0.0f, 0.0f);
		for (int i = 0; i < 3; i++) {
			Delta_b[k][i] += alpha_b[k] * (1.0f - 2.0f * pow(a_b[k][i], 2.0f) / pow(1.0f - pow(a_f[k][i], 2.0f), 2.0f));
			Delta_b[k][i] +=
				alpha_f[k] *
				(2.0f * pow(1.0f - pow(a_f[k][i], 2.0f), 2.0f) + 4.0f * a_f[k][i] * a_f[k][i] * a_b[k][i] * a_b[k][i]) /
				pow(1.0f - pow(a_f[k][i], 2.0f), 3.0f);
		}
	}
}

// The following code is only for debug purpose.
class SimpleLCG {
  public:
	SimpleLCG()
		: state(12345) {
	} // Default seed

	void seed(uint32_t seed_val) {
		// The seed cannot be 0 for this specific LCG.
		state = (seed_val == 0) ? 1 : seed_val;
	}

	// Generates the next pseudo-random integer
	uint32_t next_uint() {
		// Parameters from "Numerical Recipes"
		state = 1664525 * state + 1013904223;
		return state;
	}

	// Generates a pseudo-random float in the range [0, 1)
	float next_float() {
		return static_cast<float>(next_uint()) / static_cast<float>(std::numeric_limits<uint32_t>::max());
	}

  private:
	uint32_t state;
};

/**
 * @brief Converts spherical coordinates (theta, phi) to a direction vector in the local coordinate system used by this BSDF.
 *
 * Based on your eval function:
 * sinTheta = v.x
 * phi = atan2(v.z, v.y)
 *
 * We can derive the inverse:
 * v.x = sin(theta)
 * v.y = cos(theta) * cos(phi)
 * v.z = cos(theta) * sin(phi)
 *
 * @param theta The elevation angle, in the range [-PI/2, PI/2].
 * @param phi The azimuthal angle, in the range [0, 2*PI].
 * @return Vector3f The normalized direction vector.
 */
Vector3f anglesToVector(float theta, float phi) {
	float cosTheta = std::cos(theta);
	return Vector3f(
		std::sin(theta),
		cosTheta * std::cos(phi),
		cosTheta * std::sin(phi));
}

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
Spectrum dEonHairBSDF::integrateBSDF(const dEonHairBSDF *hairBSDF, const Vector3f &wi,
									 int integration_theta_steps, int integration_phi_steps) {

	Spectrum accumulated_value(0.0f);
	Intersection dummy_its;

	const float d_theta_o = M_PI / (float)integration_theta_steps;
	const float d_phi_o = 2.0f * M_PI / (float)integration_phi_steps;

	// --- INTEGRATION OVER THE ENTIRE SPHERE ---
	// The combination of theta_o from -90 to +90 degrees and phi_o from 0 to 360
	// degrees covers all possible outgoing directions.

	// Loop over theta_o from -PI/2 to +PI/2 (covers both hemispheres)
	for (int j = 0; j < integration_theta_steps; ++j) {
		float theta_o = -M_PI / 2.0f + (j + 0.5f) * d_theta_o;
		float cos_theta_o = std::cos(theta_o);

		// Loop over phi_o from 0 to 2*PI (covers full azimuthal rotation)
		for (int k = 0; k < integration_phi_steps; ++k) {
			float phi_o = (k + 0.5f) * d_phi_o;

			Vector3f wo = anglesToVector(theta_o, phi_o);

			BSDFSamplingRecord bRec(dummy_its, wi, wo);
			Spectrum fs_val = hairBSDF->eval(bRec, EMeasure::ESolidAngle);

			// Add the contribution of this sample to the integral.
			// The differential solid angle d(omega) is cos(theta_o) * d_theta_o * d_phi_o.
			// This term is essential for correct integration on a sphere.
			accumulated_value += fs_val * cos_theta_o;
		}
	}

	// Final multiplication by the differential area element.
	return accumulated_value * d_theta_o * d_phi_o;
}

void dEonHairBSDF::generateIntegratedHairBSDFTable(
	const dEonHairBSDF *hairBSDF,
	const std::string &filename,
	int theta_i_steps,
	int integration_theta_o_steps,
	int integration_phi_o_steps) {

	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Error: Could not open file for writing: " << filename << std::endl;
		return;
	}

	std::cout << "Starting to generate integrated table (text format)..." << std::endl;

	const float min_theta = -M_PI / 2.0f;
	const float max_theta = M_PI / 2.0f;

	// --- Write Header ---
	outFile << theta_i_steps << std::endl;
	outFile << std::fixed << std::setprecision(8) << min_theta << std::endl;
	outFile << max_theta << std::endl;

	// --- Generate and Write Data Body ---
	for (int i = 0; i < theta_i_steps; ++i) {
		float theta_i = min_theta + (float)i / (float)(theta_i_steps - 1) * (max_theta - min_theta);
		Vector3f wi = anglesToVector(theta_i, 0.0f);

		Spectrum integral_value = integrateBSDF(hairBSDF, wi, integration_theta_o_steps, integration_phi_o_steps);

		const float *val = integral_value.data();
		outFile << val[0] << " " << val[1] << " " << val[2] << std::endl;

		// Progress reporting
		if ((i + 1) % (theta_i_steps / 10) == 0 || i == theta_i_steps - 1) {
			std::cout << "\rProgress: " << (int)((float)(i + 1) / theta_i_steps * 100.0f) << "%" << std::flush;
		}
	}

	outFile.close();
	std::cout << "\nText table generation complete." << std::endl;
}

/**
 * @brief Generates a lookup table for the dEonHairBSDF::eval function and outputs it to a file.
 *
 * This function iterates through the incident angle theta_i, the outgoing angle theta_o, and the azimuthal difference angle phi_d,
 * calculates the corresponding BSDF value, and writes the results to a file in CSV format.
 *
 * @param hairBSDF A pointer to the dEonHairBSDF object to be evaluated.
 * @param filename The name of the output file.
 * @param theta_steps The number of sampling steps for the theta angle (elevation).
 * @param phi_steps The number of sampling steps for the phi angle (azimuthal difference).
 */
void generateHairBSDFTable(
	const dEonHairBSDF *hairBSDF,
	const std::string &filename,
	int theta_steps = 90, // theta from -90 degrees to +90 degrees, one point every 2 degrees
	int phi_steps = 180	  // phi_d from 0 degrees to 360 degrees, one point every 2 degrees
) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
		return;
	}

	std::cout << "Starting to generate Hair BSDF table..." << std::endl;
	std::cout << "Resolution: theta_steps=" << theta_steps << ", phi_steps=" << phi_steps << std::endl;
	std::cout << "Output file: " << filename << std::endl;

	// Write the CSV header
	outFile << "# theta_i (rad),theta_o (rad),phi_d (rad),R,G,B" << std::endl;

	// Set output precision
	outFile << std::fixed << std::setprecision(3);

	// Construct a dummy Intersection object, as it's not used in the eval function
	Intersection dummy_its;

	// Iterate over theta_i
	for (int i = 0; i < theta_steps; ++i) {
		float theta_i = -M_PI / 2.0f + (float)i / (float)(theta_steps - 1) * M_PI;

		// Print progress
		std::cout << "Progress: " << std::fixed << std::setprecision(1)
				  << (float)i / (float)theta_steps * 100.0f << "%\r" << std::flush;

		// Iterate over theta_o
		for (int j = 0; j < theta_steps; ++j) {
			float theta_o = -M_PI / 2.0f + (float)j / (float)(theta_steps - 1) * M_PI;

			// Iterate over phi_d (azimuthal difference)
			for (int k = 0; k < phi_steps; ++k) {
				float phi_d = (float)k / (float)(phi_steps - 1) * (2.0f * M_PI);

				// Exploit azimuthal symmetry, fix phi_i = 0
				float phi_i = 0.0f;
				float phi_o = phi_d;

				// Convert angles to direction vectors
				Vector3f wi = anglesToVector(theta_i, phi_i);
				Vector3f wo = anglesToVector(theta_o, phi_o);

				// Construct BSDFSamplingRecord
				// Here we use the constructor that takes wi and wo
				BSDFSamplingRecord bRec(dummy_its, wi, wo);

				// Call the eval function
				// The EMeasure parameter is not used in your eval implementation, so any value can be passed
				Spectrum result = hairBSDF->eval(bRec, EMeasure::ESolidAngle);
				if (result[0] > 1e-4) {
					// Write the result to the file
					// Assuming Spectrum can access RGB components via the [] operator
					outFile << theta_i << "," << theta_o << "," << phi_d << ","
							<< result[0] << "," << result[1] << "," << result[2] << std::endl;
				}
			}
		}
	}

	outFile.close();
	std::cout << "\nTable generation complete. Data saved to " << filename << std::endl;
}

void generateHairPDFTable(
	const dEonHairBSDF *hairBSDF,
	const std::string &filename,
	int theta_steps = 90,
	int phi_steps = 180) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
		return;
	}

	std::cout << "Starting to generate Hair PDF table..." << std::endl;
	std::cout << "Resolution: theta_steps=" << theta_steps << ", phi_steps=" << phi_steps << std::endl;
	std::cout << "Output file: " << filename << std::endl;

	// Write the CSV header
	outFile << "# theta_i (rad),theta_o (rad),phi_d (rad),PDF" << std::endl;

	// Set output precision to ensure comparable values
	outFile << std::fixed << std::setprecision(3);

	// Construct a dummy Intersection object
	Intersection dummy_its;

	// Iterate over theta_i
	for (int i = 0; i < theta_steps; ++i) {
		float theta_i = -M_PI / 2.0f + (float)i / (float)(theta_steps - 1) * M_PI;

		// Print progress
		std::cout << "Progress: " << std::fixed << std::setprecision(1)
				  << (float)i / (float)theta_steps * 100.0f << "%\r" << std::flush;

		// Iterate over theta_o
		for (int j = 0; j < theta_steps; ++j) {
			float theta_o = -M_PI / 2.0f + (float)j / (float)(theta_steps - 1) * M_PI;

			// Iterate over phi_d (azimuthal difference)
			for (int k = 0; k < phi_steps; ++k) {
				float phi_d = (float)k / (float)(phi_steps - 1) * (2.0f * M_PI);

				// Exploit azimuthal symmetry, fix phi_i = 0
				float phi_i = 0.0f;
				float phi_o = phi_d;

				// Convert angles to direction vectors
				Vector3f wi = anglesToVector(theta_i, phi_i);
				Vector3f wo = anglesToVector(theta_o, phi_o);

				// Construct BSDFSamplingRecord
				BSDFSamplingRecord bRec(dummy_its, wi, wo);

				// Call the pdf function
				// The EMeasure parameter is not used in your pdf implementation, so any value can be passed
				float pdf_value = hairBSDF->pdf(bRec, EMeasure::ESolidAngle);

				// Write the result to the file
				outFile << theta_i << "," << theta_o << "," << phi_d << "," << pdf_value << std::endl;
			}
		}
	}

	outFile.close();
	std::cout << "\nPDF table generation complete. Data saved to " << filename << std::endl;
}

void vectorToAngles_Mitsuba(const Vector3f &v, float &theta, float &phi) {
	theta = math::safe_asin(v.x);
	phi = std::atan2(v.z, v.y);
}
void drive_mitsuba_sample_test(const dEonHairBSDF *bsdf, int num_tests, const std::string filename) {
	std::ofstream outFile(filename);
	if (!outFile.is_open()) {
		std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
		return;
	}
	outFile << std::fixed << std::setprecision(3);

	outFile << "# TestID,theta_i,phi_i,xiN.x,xiN.y,xiM.x,xiM.y,wo.x,wo.y,wo.z,theta_o,phi_o" << std::endl;

	SimpleLCG lcg;
	lcg.seed(54321); // Use a fixed seed for reproducibility

	for (int i = 0; i < num_tests; ++i) {
		std::cout << "Progress: " << std::fixed << std::setprecision(1) << (float)i / (float)num_tests * 100.0f
				  << "%\r" << std::flush;

		// 1. Generate consistent test inputs using the LCG
		float u1 = lcg.next_float();
		float u2 = lcg.next_float();
		Point2f controlled_xiN(lcg.next_float(), lcg.next_float());
		Point2f controlled_xiM(lcg.next_float(), lcg.next_float());

		// Map u1, u2 to spherical coordinates for wi
		float theta_i = std::acos(1.0f - 2.0f * u1) - (M_PI / 2.0f); // Even distribution over sphere surface -> theta
		float phi_i = 2.0f * M_PI * u2;

		Vector3f wi = anglesToVector(theta_i, phi_i);

		// 2. Call the test function
		Vector3f wo = bsdf->test_sample(wi, controlled_xiN, controlled_xiM);

		// 3. Convert result back to spherical coordinates for comparison
		float theta_o, phi_o;
		vectorToAngles_Mitsuba(wo, theta_o, phi_o);

		// 4. Print results in CSV format
		outFile << i << ","
				<< theta_i << "," << phi_i << ","
				<< controlled_xiN.x << "," << controlled_xiN.y << ","
				<< controlled_xiM.x << "," << controlled_xiM.y << ","
				<< theta_o << "," << phi_o << std::endl;
	}
}

Vector3f dEonHairBSDF::test_sample(
	const Vector3f &wi,
	const Point2f &controlled_xiN,
	const Point2f &controlled_xiM) const {
	// --- Our controlled inputs ---
	const Point2 &xiN = controlled_xiN;
	const Point2 &xiM = controlled_xiM;
	// ---

	Float sinThetaI = wi.x;
	Float cosThetaI = math::trigInverse(sinThetaI);
	Float thetaI = math::safe_asin(sinThetaI);
	Float phiI = std::atan2(wi.z, wi.y);

	Float thetaIR = thetaI - 2.0f * _scaleAngleRad;
	Float thetaITT = thetaI + _scaleAngleRad;
	Float thetaITRT = thetaI + 4.0f * _scaleAngleRad;

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
	Float sinPhiO = std::sin(phi + phiI);
	Float cosPhiO = std::cos(phi + phiI);

	Vector3 wo = Vector3(
		sinThetaO,
		cosPhiO * cosThetaO,
		sinPhiO * cosThetaO);

	return wo;
}

MTS_IMPLEMENT_CLASS_API(dEonHairBSDF, false, BSDF, DEON_API)
MTS_EXPORT_PLUGIN(dEonHairBSDF, "Marschner Hair Shading Model");
MTS_NAMESPACE_END
