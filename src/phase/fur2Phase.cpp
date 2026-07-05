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
// #pragma optimize("", off)
#include "fur2Phase.h"
#include <map>
#include <mitsuba/core/warp.h>
#include <mitsuba/hw/basicshader.h>
#include <vector>

using namespace std;

MTS_NAMESPACE_BEGIN

void printTablesToFiles(
	const std::vector<Vector3> &full_table, const std::string &full_filename) {
	const size_t theta_i_steps = full_table.size();
	if (theta_i_steps == 0) {
		std::cout << "Warning: table is empty; nothing to write." << std::endl;
		return;
	}

	const Float min_theta_i_rad = -M_PI / 2.0f;
	const Float max_theta_i_rad = M_PI / 2.0f;
	const Float rad_to_deg = 180.0f / M_PI;

	std::ofstream full_file(full_filename);
	if (!full_file.is_open()) {
		std::cerr << "Error: could not open file " << full_filename << std::endl;
		return;
	}

	std::cout << "Writing " << full_filename << "..." << std::endl;
	for (size_t i = 0; i < theta_i_steps; ++i) {
		Float theta_i_rad = min_theta_i_rad + (Float)i / (Float)(theta_i_steps - 1) * (max_theta_i_rad - min_theta_i_rad);
		Float theta_i_deg = theta_i_rad * rad_to_deg;
		const Vector3 &vec = full_table[i];

		full_file << theta_i_deg << " " << vec.x << " " << vec.y << " " << vec.z << "\n";
	}
	std::cout << "Done." << std::endl;
}

Fur2Phase::Fur2Phase(const Properties &props)
	: SpectrumPhaseFunction(props) {
	m_totalScatteredEnergy = std::make_shared<std::vector<Vector3>>();

	// If color texture is given, read and convert it into absorption on the fly
	colorTexture = new ConstantSpectrumTexture(props.getSpectrum("color", Spectrum(1.0f)));

	useColorTexture = props.getBoolean("useColorTexture", false);
	medulla_ratio = props.getFloat("medulla_ratio", 0.0);
	eta = props.getFloat("eta", 1.55);
	beta_m = props.getFloat("beta_m", 0.055);
	beta_n = props.getFloat("beta_n", 0.055);
	alpha = props.getFloat("alpha", 3.0);
	g = props.getFloat("g", 0.0);
	layers = props.getFloat("layers", 0.5);
	absorptionFactor = props.getFloat("absorptionFactor", 4.0);
	absorption_outer = props.getVector("absorption_outer", Vector3(0.0f, 0.0f, 0.0f));
	absorption_inner = props.getFloat("absorption_inner", 0.0f);
	scattering_inner = props.getFloat("scattering_inner", 0.5f);
	furMaskCode = props.getInteger("furMaskCode", 31);

	// Render mode:
	//     "near":      near field only
	//     "far":       far field only
	//     "nearfar":   near field for first bounce, far field for others
	renderMode = props.getString("renderMode", "far");
	if (renderMode != "near" && renderMode != "far" && renderMode != "nearfar" && renderMode != "midfar")
		printf("Invalid input for render mode!\n");

	// Color scale:
	//     > 0:             scale by colorScale
	// Color gamma:
	//     > 0:             gamma by colorScale
	colorScale = props.getFloat("colorScale", 2.0);
	colorGamma = props.getFloat("colorGamma", 1.0);

	// Tilt should always be negative for physical correctness.
	alpha = -degToRad(alpha);

	initialize("medulla");

	generateIntegratedHairBSDFTable();
	printTablesToFiles(*m_totalScatteredEnergy, "fur2PhaseTable.txt");
}

Fur2Phase::Fur2Phase(Stream *stream, InstanceManager *manager)
	: SpectrumPhaseFunction(stream, manager) {
	colorTexture = static_cast<Texture *>(manager->getInstance(stream));

	//        perspectiveCamera = static_cast<PerspectiveCamera *>(manager->getInstance(stream));

	useColorTexture = stream->readBool();
	medulla_ratio = stream->readFloat();
	eta = stream->readFloat();
	beta_m = stream->readFloat();
	beta_n = stream->readFloat();
	alpha = stream->readFloat();
	g = stream->readFloat();
	layers = stream->readFloat();
	absorptionFactor = stream->readFloat();
	absorption_outer.x = stream->readFloat();
	absorption_outer.y = stream->readFloat();
	absorption_outer.z = stream->readFloat();
	absorption_inner = stream->readFloat();
	scattering_inner = stream->readFloat();

	renderMode = stream->readString();
	colorScale = stream->readFloat();
	colorGamma = stream->readFloat();

	initialize("medulla");

	configure();
}

void Fur2Phase::configure() {
	/* Verify the input parameter and fix them if necessary */
	// colorTexture = ensureEnergyConservation(colorTexture, "color", 1.0f);

	// m_components.clear();
	// m_components.push_back(EFUR2 | EUsesSampler);
	// if (colorTexture->getMaximum().max() > 0)
	//	m_components.push_back(EDiffuseReflection | EFrontSide | (colorTexture->isConstant() ? 0 : ESpatiallyVarying));
	// m_usesRayDifferentials = colorTexture->usesRayDifferentials();

	// BSDF::configure();
}

Spectrum Fur2Phase::getDiffuseReflectance(const Intersection &its) const {
	return colorTexture->eval(its);
}

Vector3 Fur2Phase::normalizeSafe(Vector3 v) const {
	if (v.length() == 0.0)
		return v;
	else
		return normalize(v);
}

Spectrum Fur2Phase::eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
	Vector3 absorption = absorption_outer;
	if (useColorTexture && medulla_ratio < 1) {
		Spectrum color = colorTexture->eval(bRec.its, false);
		Float color_r, color_g, color_b;
		color.toLinearRGB(color_r, color_g, color_b);

		color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
		color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
		color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

		Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
		Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
		Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
		Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

		absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
		absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
		absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

		absorption = Vector3(absorption_r, absorption_g, absorption_b);
	}

	Vector3 wi = bRec.its.toWorld(bRec.wi);
	Vector3 wo = bRec.its.toWorld(bRec.wo);
	const Vector3 &u = bRec.its.geoFrame.s;
	const Vector3 &v = bRec.its.geoFrame.t;
	const Vector3 &w = bRec.its.geoFrame.n;

	Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
	Float theta_r = M_PI / 2.0f - acosSafe(dot(u, wo));
	theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
	theta_r = clamp(theta_r, degToRad(-89.5), degToRad(89.5));

	Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
	Float phi_r = acosSafe(dot(normalizeSafe(wo - dot(u, wo) * u), v)) * (dot(w, wo) > 0 ? 1 : -1);
	Float phi = phi_r - phi_i;
	regularizePhi(phi);

	FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
						medulla_ratio, absorption_inner, scattering_inner, g, layers, bRec.sampler);

	Vector3 R, TT, TRT, TTs, TRTs;

	Float ms = furSolver.Ms(theta_i - alpha, theta_r - alpha, phi);

	if (renderMode == "near" || (renderMode == "nearfar" && bRec.depth == 1)) {
		Log(EError, "Not supported!");
		Float h = bRec.its.storage[1];
		R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.Nu(0, h, phi, theta_i);
		TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.Nu(1, h, phi, theta_i);
		TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.Nu(2, h, phi, theta_i);
		TTs = ms * furSolver.Ns(1, h, phi, theta_i);
		TRTs = ms * furSolver.Ns(2, h, phi, theta_i);
	} else if (renderMode == "midfar" && bRec.depth == 1) {
		Log(EError, "Not supported!");
		Point origin = perspectiveCamera->s_origin;
		Vector3 lookAt = perspectiveCamera->s_lookAt;

		Point p = bRec.its.p;
		Vector3 pVec = p - origin;
		Float projDist = dot(pVec, lookAt);
		Float fovX = perspectiveCamera->getXFov() / 180.0 * M_PI;
		Float covX = tan(fovX / 2.0) * projDist * 2.0 * PIXEL_COVERAGE_SCALE;
		int resX = perspectiveCamera->getFilm()->getSize().x;
		Float pixCov = covX / resX;
		Float furRadius = bRec.its.storage[2];
		Vector3 pVecNorm = normalize(pVec);
		Float pixCovProj = pixCov * dot(pVecNorm, lookAt);
		Float hRange = pixCovProj / furRadius;

		Float h = bRec.its.storage[1];
		Float hMin = clamp(h - hRange / 2.0, -1.0 + EPS, 1.0 - EPS);
		Float hMax = clamp(h + hRange / 2.0, -1.0 + EPS, 1.0 - EPS);

		//            if (hMin != -0.999 || hMax != 0.999)
		//                printf("%f %f\n", hMin, hMax);

		if (hMin > 0.0) {
			R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.NuMid(0, hMin, hMax, phi, theta_i);
			TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.NuMid(1, hMin, hMax, phi, theta_i);
			TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.NuMid(2, hMin, hMax, phi, theta_i);
			TTs = ms * max(discardLarger(furSolver.NsMid(1, hMin, hMax, phi, theta_i), 10.0f), 0.0f);
			TRTs = ms * max(discardLarger(furSolver.NsMid(2, hMin, hMax, phi, theta_i), 10.0f), 0.0f);
		} else if (hMax < 0.0) {
			R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.NuMid(0, -hMax, -hMin, -phi, theta_i);
			TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.NuMid(1, -hMax, -hMin, -phi, theta_i);
			TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.NuMid(2, -hMax, -hMin, -phi, theta_i);
			TTs = ms * max(discardLarger(furSolver.NsMid(1, -hMax, -hMin, -phi, theta_i), 10.0f), 0.0f);
			TRTs = ms * max(discardLarger(furSolver.NsMid(2, -hMax, -hMin, -phi, theta_i), 10.0f), 0.0f);
		} else {
			R = furSolver.Mu(0, theta_i, theta_r, phi) * (furSolver.NuMid(0, 0, hMax, phi, theta_i) * hMax / (hMax - hMin) + furSolver.NuMid(0, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin));
			TT = furSolver.Mu(1, theta_i, theta_r, phi) * (furSolver.NuMid(1, 0, hMax, phi, theta_i) * hMax / (hMax - hMin) + furSolver.NuMid(1, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin));
			TRT = furSolver.Mu(2, theta_i, theta_r, phi) * (furSolver.NuMid(2, 0, hMax, phi, theta_i) * hMax / (hMax - hMin) + furSolver.NuMid(2, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin));
			TTs = ms * (max(discardLarger(furSolver.NsMid(1, 0, hMax, phi, theta_i) * hMax / (hMax - hMin), 10.0f), 0.0f) + max(discardLarger(furSolver.NsMid(1, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin), 10.0f), 0.0f));
			TRTs = ms * (max(discardLarger(furSolver.NsMid(2, 0, hMax, phi, theta_i) * hMax / (hMax - hMin), 10.0f), 0.0f) + max(discardLarger(furSolver.NsMid(2, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin), 10.0f), 0.0f));
		}
		////printf("%f %f %f\n", R.x, R.y, R.z);
		//            TT = Vector3(0.0, 0.0, 0.0);
		//            TRT = Vector3(0.0, 0.0, 0.0);
		//            TTs = Vector3(0.0, 0.0, 0.0);
		//            TRTs = Vector3(0.0, 0.0, 0.0);
	} else {
		R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.Nu(0, phi, theta_i);
		TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.Nu(1, phi, theta_i);
		TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.Nu(2, phi, theta_i);
		TTs = ms * max(discardLarger(furSolver.Ns(1, phi, theta_i), 10.0f), 0.0f);
		TRTs = ms * max(discardLarger(furSolver.Ns(2, phi, theta_i), 10.0f), 0.0f);

		//            if (max(discardLarger(furSolver.Ns(1, phi, theta_i), 10.0f), 0.0f).length() > 1.0 ||
		//                max(discardLarger(furSolver.Ns(2, phi, theta_i), 10.0f), 0.0f).length() > 1.0)
		//            printf("R = %f %f %f\nTT=%f %f %f\nTRT=%f %f %f\nTTs=%f %f %f\nTRTs=%f %f %f\n\n",
		//                   R.x, R.y, R.z,
		//                   TT.x, TT.y, TT.z,
		//                   TRT.x, TRT.y, TRT.z,
		//                   TTs.x, TTs.y, TTs.z,
		//                   TRTs.x, TRTs.y, TRTs.z);
	}
	R = furMaskCode & 1 ? R : Vector3(Float(0));
	TT = (furMaskCode >> 1) & 1 ? TT : Vector3(Float(0));
	TRT = (furMaskCode >> 2) & 1 ? TRT : Vector3(Float(0));
	TTs = (furMaskCode >> 3) & 1 ? TTs : Vector3(Float(0));
	TRTs = (furMaskCode >> 4) & 1 ? TRTs : Vector3(Float(0));

	Vector3 bsdfCosine;
	if (bRec.component == -1)
		bsdfCosine = (R + TT + TRT + TTs + TRTs) / cos(theta_r);
	else if (bRec.component == 0)
		bsdfCosine = (R + TT + TRT) / cos(theta_r);
	else
		bsdfCosine = (TTs + TRTs) / cos(theta_r);

	if (std::isnan(bsdfCosine[0]) || std::isnan(bsdfCosine[1]) || std::isnan(bsdfCosine[2]))
		return Spectrum(0.0f);

	if (std::isinf(bsdfCosine[0]) || std::isinf(bsdfCosine[1]) || std::isinf(bsdfCosine[2]))
		return Spectrum(0.0f);

	for (int i = 0; i < 3; i++)
		bsdfCosine[i] = std::min(std::max(Float(0.0f), bsdfCosine[i]), Float(1e3f));

	Spectrum retVal(0.0f);
	retVal.fromLinearRGB(bsdfCosine[0], bsdfCosine[1], bsdfCosine[2]);

	// Vector3 ret = Vector3(1, 1, 1)*1e5; retVal.fromLinearRGB(ret[0], ret[1], ret[2]); // for debug

	return retVal;
}

Float Fur2Phase::normrnd(const BSDFSamplingRecord &bRec, Float mu, Float sigma) const {
	double U = bRec.sampler->next1D();
	double V = bRec.sampler->next1D();
	double X = sqrt(-2.0 * log(U)) * cos(2.0 * M_PI * V);
	return X * sigma + mu;
}

Float Fur2Phase::pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
	Vector3 absorption = absorption_outer;
	if (useColorTexture && medulla_ratio < 1) {
		Spectrum color = colorTexture->eval(bRec.its, false);
		Float color_r, color_g, color_b;
		color.toLinearRGB(color_r, color_g, color_b);

		color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
		color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
		color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

		Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
		Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
		Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
		Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

		absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
		absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
		absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

		absorption = Vector3(absorption_r, absorption_g, absorption_b);
	}

	Vector3 wi = bRec.its.toWorld(bRec.wi);
	Vector3 wo = bRec.its.toWorld(bRec.wo);
	const Vector3 &u = bRec.its.geoFrame.s;
	const Vector3 &v = bRec.its.geoFrame.t;
	const Vector3 &w = bRec.its.geoFrame.n;

	Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
	Float theta_r = M_PI / 2.0f - acosSafe(dot(u, wo));
	theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
	theta_r = clamp(theta_r, degToRad(-89.5), degToRad(89.5));

	Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
	Float phi_r = acosSafe(dot(normalizeSafe(wo - dot(u, wo) * u), v)) * (dot(w, wo) > 0 ? 1 : -1);
	Float phi = phi_r - phi_i;
	regularizePhi(phi);

	FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
						medulla_ratio, absorption_inner, scattering_inner, g, layers, bRec.sampler);

	Float h;
	if (renderMode == "near" || (renderMode == "nearfar" && bRec.depth == 1)) {
		h = bRec.its.storage[1];
	} else if (renderMode == "midfar") {
		Point origin = perspectiveCamera->s_origin;
		Vector3 lookAt = perspectiveCamera->s_lookAt;
		Point p = bRec.its.p;
		Vector3 pVec = p - origin;
		Float projDist = dot(pVec, lookAt);
		Float fovX = perspectiveCamera->getXFov() / 180.0 * M_PI;
		Float covX = tan(fovX / 2.0) * projDist * 2.0 * PIXEL_COVERAGE_SCALE;
		int resX = perspectiveCamera->getFilm()->getSize().x;
		Float pixCov = covX / resX;
		Float furRadius = bRec.its.storage[2];
		Vector3 pVecNorm = normalize(pVec);
		Float pixCovProj = pixCov * dot(pVecNorm, lookAt);
		Float hRange = pixCovProj / furRadius;

		Float hMid = bRec.its.storage[1];
		Float hMin = clamp(hMid - hRange / 2.0, -1.0 + EPS, 1.0 - EPS);
		Float hMax = clamp(hMid + hRange / 2.0, -1.0 + EPS, 1.0 - EPS);

		h = bRec.sampler->next1D() * (hMax - hMin) + hMin;
	} else {
		h = (bRec.sampler->next1D() * 2.0 - 1.0) * (1.0 - EPS);
	}

	return std::max(Float(0.0f), furSolver.pdfMN(h, theta_i, theta_r, phi));
}

Spectrum Fur2Phase::sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
	// bRec.eta = Float(1.f);
	Vector3 absorption = absorption_outer;
	if (useColorTexture && medulla_ratio < 1) {
		Spectrum color = colorTexture->eval(bRec.its, false);
		Float color_r, color_g, color_b;
		color.toLinearRGB(color_r, color_g, color_b);

		color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
		color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
		color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

		Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
		Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
		Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
		Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

		absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
		absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
		absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

		absorption = Vector3(absorption_r, absorption_g, absorption_b);
	}

	bRec.sampledComponent = 0;
	bRec.sampledType = BSDF::EDiffuseReflection;

	Vector3 wi = bRec.its.toWorld(bRec.wi);
	const Vector3 &u = bRec.its.geoFrame.s;
	const Vector3 &v = bRec.its.geoFrame.t;
	const Vector3 &w = bRec.its.geoFrame.n;

	Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
	theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));

	Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);

	FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
						medulla_ratio, absorption_inner, scattering_inner, g, layers, bRec.sampler);

	// Choose a lobe to sample.
	Float h;
	if (renderMode == "near" || (renderMode == "nearfar" && bRec.depth == 1)) {
		h = bRec.its.storage[1];
	} else if (renderMode == "midfar") {
		Point origin = perspectiveCamera->s_origin;
		Vector3 lookAt = perspectiveCamera->s_lookAt;
		Point p = bRec.its.p;
		Vector3 pVec = p - origin;
		Float projDist = dot(pVec, lookAt);
		Float fovX = perspectiveCamera->getXFov() / 180.0 * M_PI;
		Float covX = tan(fovX / 2.0) * projDist * 2.0 * PIXEL_COVERAGE_SCALE;
		int resX = perspectiveCamera->getFilm()->getSize().x;
		Float pixCov = covX / resX;
		Float furRadius = bRec.its.storage[2];
		Vector3 pVecNorm = normalize(pVec);
		Float pixCovProj = pixCov * dot(pVecNorm, lookAt);
		Float hRange = pixCovProj / furRadius;

		Float hMid = bRec.its.storage[1];
		Float hMin = clamp(hMid - hRange / 2.0, -1.0 + EPS, 1.0 - EPS);
		Float hMax = clamp(hMid + hRange / 2.0, -1.0 + EPS, 1.0 - EPS);

		h = bRec.sampler->next1D() * (hMax - hMin) + hMin;
	} else {
		h = (bRec.sampler->next1D() * 2.0 - 1.0) * (1.0 - EPS);
	}

	int lobe;
	Float pdfLobe;
	furSolver.chooseLobe(theta_i, h, lobe, pdfLobe);

	bRec.sampledComponent = lobe < 3 ? 0 : 1;

	if (lobe == -1) {
		pdf = 0.0f;
		return Spectrum(0.0f);
	}

	if (((furMaskCode >> lobe) & 1) == 0) {
		pdf = 0.0f;
		return Spectrum(0.0f);
	}

	// Sample the chosen lobe.
	Float theta_o, phi_o;
	Float pdfM, pdfN;
	Float weightM;
	Vector3 weightN;
	furSolver.sampleLobe(lobe, h, theta_i, phi_i, theta_o, pdfM, weightM, phi_o, pdfN, weightN);

	regularizePhi(phi_o);
	regularizeTheta(theta_o);
	theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));

	// Convert angles to outgoing direction.
	Float ul = sin(theta_o);
	Float vl = cos(theta_o) * cos(phi_o);
	Float wl = cos(theta_o) * sin(phi_o);
	bRec.wo = bRec.its.toLocal(u * ul + v * vl + w * wl);

	// Calculate the final sampling weight and pdf.
	Vector3 weight = weightM * weightN / pdfLobe;

	if (std::isnan(weight[0]) || std::isnan(weight[1]) || std::isnan(weight[2]))
		return Spectrum(0.0f);

	if (std::isinf(weight[0]) || std::isinf(weight[1]) || std::isinf(weight[2]))
		return Spectrum(0.0f);

	weight[0] = std::max(Float(0.0f), std::min(weight[0], Float(2.0f)));
	weight[1] = std::max(Float(0.0f), std::min(weight[1], Float(2.0f)));
	weight[2] = std::max(Float(0.0f), std::min(weight[2], Float(2.0f)));

	// Compute the pdf. This can be more efficient.
	Float phi = phi_o - phi_i;
	pdf = pdfM * pdfN * pdfLobe / cos(theta_o);

	Spectrum retVal(0.f);
	retVal.fromLinearRGB(weight[0], weight[1], weight[2]);
	return retVal;
}

Spectrum Fur2Phase::sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
	Float pdf;
	return Fur2Phase::sample(bRec, pdf, sample);
}

PFPrecalcSample Fur2Phase::sampleWithLocalWiForIntegration(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const {
	const Vector3 u(1.0f, 0.0f, 0.0f);
	const Vector3 v(0.0f, 1.0f, 0.0f);
	const Vector3 w(0.0f, 0.0f, 1.0f);

	Float theta_i = M_PI / 2.0f - acosSafe(wi.x);
	theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));

	// Project wi onto the YZ plane and compute its azimuth.
	Vector3 wi_perp = wi - wi.x * u; // (0, wi.y, wi.z)
	Vector3 wi_perp_norm = normalizeSafe(wi_perp);

	// dot(wi_perp_norm, v) -> wi_perp_norm.y
	// dot(w, wi) > 0 -> wi.z > 0
	Float phi_i = acosSafe(wi_perp_norm.y) * (wi.z > 0 ? 1 : -1);

	Vector3 absorption = absorption_outer;

	FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
						medulla_ratio, absorption_inner, scattering_inner, g, layers, sampler);

	// Choose a lobe to sample.
	Float h;
	if (renderMode == "near" || (renderMode == "nearfar" /*&& bRec.depth == 1*/)) {
		Log(EError, "near and nearfar not supported");
	} else if (renderMode == "midfar") {
		Log(EError, "midfar not supported");

	} else {
		h = (sampler->next1D() * 2.0 - 1.0) * (1.0 - EPS);
	}

	int lobe;
	Float pdfLobe;
	furSolver.chooseLobe(theta_i, h, lobe, pdfLobe);

	PFPrecalcSample result;

	if (lobe == -1 || ((furMaskCode >> lobe) & 1) == 0) {
		pdf = 0.0f;
		wo = Vector3(0.0f);
		result.m_totalScatteredEnergy = Vector3(0.f);
		return result;
	}

	Float theta_o, phi_o;
	Float pdfM, pdfN;
	Float weightM;
	Vector3 weightN;
	furSolver.sampleLobe(lobe, h, theta_i, phi_i, theta_o, pdfM, weightM, phi_o, pdfN, weightN);

	regularizePhi(phi_o);
	regularizeTheta(theta_o);
	theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));

	// Convert sampled angles to a local outgoing direction.
	Float ul = sin(theta_o);
	Float vl = cos(theta_o) * cos(phi_o);
	Float wl = cos(theta_o) * sin(phi_o);

	wo = Vector3(ul, vl, wl);

	Vector3 weight = weightM * weightN / pdfLobe;

	if (std::isnan(weight[0]) || std::isnan(weight[1]) || std::isnan(weight[2]) ||
		std::isinf(weight[0]) || std::isinf(weight[1]) || std::isinf(weight[2])) {
		result.m_totalScatteredEnergy = Vector3(0.f);
		return result;
	}

	weight[0] = std::max(Float(0.0f), std::min(weight[0], Float(2.0f)));
	weight[1] = std::max(Float(0.0f), std::min(weight[1], Float(2.0f)));
	weight[2] = std::max(Float(0.0f), std::min(weight[2], Float(2.0f)));

	// Convert from theta/phi density to solid-angle density.
	pdf = pdfM * pdfN * pdfLobe / cos(theta_o);

	Spectrum retVal(0.f);
	retVal.fromLinearRGB(weight[0], weight[1], weight[2]);

	result.m_totalScatteredEnergy = weight;

	return result;
}

Float Fur2Phase::eval(const PhaseFunctionSamplingRecord &pRec) const {
	Log(EError, "Not supported.");
}

Float Fur2Phase::sample(PhaseFunctionSamplingRecord &pRec,
						Sampler *sampler) const {
	Log(EError, "Not supported.");
}

Float Fur2Phase::sample(PhaseFunctionSamplingRecord &pRec,
						Float &pdf, Sampler *sampler) const {
	Log(EError, "Not supported.");
}

Vector3 Fur2Phase::evalRGB(const PhaseFunctionSamplingRecord &pRec, const Float sample1D) const {
	if (pRec.mRec.orientation.isZero()) {
		return Vector3(0.0f);
	}

	Frame frame = Frame::fromS(pRec.mRec.orientation);
	Vector wi = frame.toLocal(pRec.wi).normalized();
	Vector wo = frame.toLocal(pRec.wo).normalized();

	Intersection its;
	its.shFrame = frame;
	its.geoFrame = frame;

	BSDFSamplingRecord bRec(its, wi, wo);

	Spectrum evalResult = eval(bRec, EMeasure::ESolidAngle);

	// Construct thetaI
	Vector3 worldWi = bRec.its.toWorld(bRec.wi);
	Vector3 worldWo = bRec.its.toWorld(bRec.wo);
	const Vector3 &u = bRec.its.geoFrame.s;
	const Vector3 &v = bRec.its.geoFrame.t;
	const Vector3 &w = bRec.its.geoFrame.n;
	Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
	Vector3 normalizationFactor = lookupTotalScatteredEnergy(theta_i);

	Vector3 ret = evalResult.toLinearRGB();
	ret = CwiseDiv(ret, normalizationFactor);

	return ret;
}

void Fur2Phase::addChild(const std::string &name, ConfigurableObject *child) {
	if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "color") {
		colorTexture = static_cast<Texture *>(child);
	} else if (child->getClass()->derivesFrom(MTS_CLASS(Sensor))) {
		//            printf("Yes\n");
		perspectiveCamera = static_cast<PerspectiveCamera *>(child);
	} else {
		SpectrumPhaseFunction::addChild(name, child);
	}
}

bool Fur2Phase::needsDirectionallyVaryingCoefficients() const {
	return true;
}

Float Fur2Phase::sigmaDir(Float cosTheta) const {
	Float sinTheta = math::trigInverse(cosTheta);
	return sinTheta;
}

Float Fur2Phase::sigmaDirMax() const {
	return 1.f;
}

Float Fur2Phase::getMeanCosine() const {
	return 0;
}

// #pragma optimize("", off)
Spectrum Fur2Phase::albedoDir(const Vector &d, const Vector &orientation) const {
	Frame frame = Frame::fromS(orientation);
	Vector wi = frame.toLocal(-d);
	Float sinThetaI = wi.x;
	Float thetaI = math::safe_asin(sinThetaI);
	Vector3 scattered = lookupTotalScatteredEnergy(thetaI);
	Spectrum result;

	// Vector newOrientiation = orientation;
	////newOrientiation.z = -newOrientiation.z;
	// Vector3 r = (newOrientiation + Vector(1.f)) / 2.f;
	// result.fromLinearRGB(r.x, r.y, r.z);
	// return result;

	result.fromLinearRGB(scattered.x, scattered.y, scattered.z);
	return result;
}
// #pragma optimize("", on)

bool Fur2Phase::needsDirectionallyVaryingAlbedo() const {
	return true;
}

Vector3 Fur2Phase::lookupTotalScatteredEnergy(Float thetaI) const {
	// --- Safety Check ---
	// Ensure the table has been generated before trying to access it.
	if (m_totalScatteredEnergy->empty()) {
		Log(EWarn, "lookupNormalizationFactor called, but the integration table is empty. Returning 0.");
		return Vector3(0.0f);
	}

	// --- Define the Domain ---
	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;
	const int table_size = m_totalScatteredEnergy->size();

	// --- Normalize and Clamp the Input ---
	// Map the input theta_i to a normalized coordinate 'u' in the range [0, 1].
	Float u = (thetaI - min_theta_i) / (max_theta_i - min_theta_i);
	u = std::max(0.0f, std::min(1.0f, u)); // Clamp to handle out-of-range inputs

	// --- Calculate Indices for Interpolation ---
	// Find the floating-point position in the table.
	// The number of intervals is (table_size - 1).
	const Float pos = u * (table_size - 1);

	// Find the two integer indices that bracket our position.
	const int idx0 = static_cast<int>(std::floor(pos));
	const int idx1 = std::min(idx0 + 1, table_size - 1); // Clamp to avoid going out of bounds

	// --- Perform Linear Interpolation ---
	// Calculate the interpolation weight 't' (the fractional part of 'pos').
	const Float t = pos - idx0;

	// Retrieve the values at the two indices.
	const Vector3 val0 = (*m_totalScatteredEnergy)[idx0];
	const Vector3 val1 = (*m_totalScatteredEnergy)[idx1];

	// Linearly interpolate between the two values and return the result.
	// Formula: (1 - t) * val0 + t * val1
	Vector3 scattered = val0 + t * (val1 - val0);
	scattered = scattered.saturate();
	Vector3 maskedScattered(0.f);
	// if (maskAlbedo) {
	//	maskedScattered = CwiseMul(scattered, rgbMask);
	// } else {
	//	maskedScattered = scattered;
	// }
	maskedScattered = scattered;
	return maskedScattered;
}

void Fur2Phase::serialize(Stream *stream, InstanceManager *manager) const {
	SpectrumPhaseFunction::serialize(stream, manager);

	manager->serialize(stream, colorTexture.get());
	//        perspectiveCamera->serialize(stream, manager);

	stream->writeBool(useColorTexture);
	stream->writeFloat(medulla_ratio);
	stream->writeFloat(eta);
	stream->writeFloat(beta_m);
	stream->writeFloat(beta_n);
	stream->writeFloat(alpha);
	stream->writeFloat(g);
	stream->writeFloat(layers);
	stream->writeFloat(absorptionFactor);
	stream->writeFloat(absorption_outer.x);
	stream->writeFloat(absorption_outer.y);
	stream->writeFloat(absorption_outer.z);
	stream->writeFloat(absorption_inner);
	stream->writeFloat(scattering_inner);

	stream->writeString(renderMode);
	stream->writeFloat(colorScale);
	stream->writeFloat(colorGamma);
}

Float Fur2Phase::getRoughness(const Intersection &its, int component) const {
	return std::numeric_limits<Float>::infinity();
}

std::string Fur2Phase::toString() const {
	std::ostringstream oss;
	oss << "SmoothDiffuse[" << endl
		<< "  id = \"" << getID() << "\"," << endl
		<< "  reflectance = " << indent(colorTexture->toString()) << endl
		<< "]";
	return oss.str();
}

PFPrecalcSample Fur2Phase::integrateBSDF(const Vector3f &wi, Sampler *sampler) const {
	PFPrecalcSample accumulated_value;

	if (sampler == nullptr) {
		Log(EError, "Sampler is null in Fur2Phase::integrateBSDF.");
		return accumulated_value;
	}

	// Reset the sampler so each integration point uses a consistent low-discrepancy sequence.
	sampler->generate(Point2i(-1));

	for (int i = 1; i <= m_intergration_sample_cnt; i++) {
		Vector3 wo;
		Float pdf;

		// Use the thread-local sampler passed by the caller.
		auto sampledValue = sampleWithLocalWiForIntegration(wi, wo, pdf, sampler);

		/* Reference code for uniform-sphere integration.
		wo = warp::squareToUniformSphere(sampler->next2D());
		auto sampledValue = evalBSDF(wi, wo);
		sampledValue = sampledValue / warp::squareToUniformSpherePdf();*/

		if (!std::isfinite(sampledValue.m_normalizationFactor)) {
			Log(EWarn, "Non-finite value detected in Fur2Phase integration.");
		}

		sampler->advance();
		accumulated_value += sampledValue;
	}

	return accumulated_value / (Float)m_intergration_sample_cnt;
}

// -----------------------------------------------------------------------
// Integrated BSDF table generation.
// Parallelized with OpenMP.
// -----------------------------------------------------------------------
void Fur2Phase::generateIntegratedHairBSDFTable() {
	Log(EInfo, "Starting to generate integrated BSDF table (Fur2Phase)...");

	if (m_theta_i_steps <= 1) {
		Log(EError, "m_theta_i_steps must be greater than 1, but is %d.", m_theta_i_steps);
		return;
	}

	// Allocate output tables.
	m_normalizationFactors.clear();
	m_normalizationFactors.resize(m_theta_i_steps);
	m_totalScatteredEnergy->clear();
	m_totalScatteredEnergy->resize(m_theta_i_steps);

	const Float min_theta_i = -M_PI / 2.0f;
	const Float max_theta_i = M_PI / 2.0f;

	// Create the prototype sampler once on the main thread.
	auto pluginManager = PluginManager::getInstance();
	Properties samplerProperties("ldsampler");
	samplerProperties.setInteger("sampleCount", this->m_intergration_sample_cnt);

	ref<ConfigurableObject> protoObj = pluginManager->createObject(MTS_CLASS(Sampler), samplerProperties);
	ref<Sampler> prototypeSampler = dynamic_cast<Sampler *>(protoObj.get());

	if (!prototypeSampler) {
		Log(EError, "Failed to create prototype sampler.");
		return;
	}

	// Atomic counter for progress logging.
	std::atomic<int> processed_count(0);

// OpenMP parallel region.
#pragma omp parallel
	{
		// Each thread receives its own sampler clone to avoid shared sampler state.
		ref<Sampler> threadSampler = dynamic_cast<Sampler *>(prototypeSampler->clone().get());

		if (!threadSampler) {
#pragma omp critical
			Log(EError, "Failed to clone sampler in thread %d", omp_get_thread_num());
		}

		// Parallel integration over incident angles.
		// Dynamic scheduling balances angles whose integrations take different times.
#pragma omp for schedule(dynamic)
		for (int i = 0; i < m_theta_i_steps; ++i) {
			Float theta_i = min_theta_i + (Float)i / (Float)(m_theta_i_steps - 1) * (max_theta_i - min_theta_i);

			// The azimuth is fixed to zero due to azimuthal symmetry.
			Vector3f wi = BCSDFAngleToVector(theta_i, 0.0f);

			// Integrate with the thread-local sampler.
			auto integral_value = integrateBSDF(wi, threadSampler);

			// Each loop iteration writes a unique index.
			m_normalizationFactors[i] = integral_value.m_normalizationFactor;
			(*m_totalScatteredEnergy)[i] = integral_value.m_totalScatteredEnergy;

			// Progress logging.
			int current = ++processed_count;
			// Log every 10% of completed work.
			int step = m_theta_i_steps / 10;
			if (step > 0 && (current % step == 0 || current == m_theta_i_steps)) {
#pragma omp critical
				{
					Log(EInfo, "Progress: %d%%", (int)((Float)current / m_theta_i_steps * 100.0f));
				}
			}
		}
		// threadSampler is released when it leaves scope.
	}

	Log(EInfo, "Integrated BSDF table generation complete. Table size: %zu.", m_normalizationFactors.size());
}

// ================ Hardware shader implementation ================

class Fur2Shader : public Shader {
  public:
	Fur2Shader(Renderer *renderer, const Texture *reflectance)
		: Shader(renderer, EBSDFShader), colorTexture(reflectance) {
		colorTextureShader = renderer->registerShaderForResource(colorTexture.get());
	}

	bool isComplete() const {
		return colorTextureShader.get() != NULL;
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(colorTexture.get());
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(colorTextureShader.get());
	}

	void generateCode(std::ostringstream &oss,
					  const std::string &evalName,
					  const std::vector<std::string> &depNames) const {
		oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "        return vec3(0.0);" << endl
			<< "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    return " << evalName << "(uv, wi, wo);" << endl
			<< "}" << endl;
	}

	MTS_DECLARE_CLASS()
  private:
	ref<const Texture> colorTexture;
	ref<Shader> colorTextureShader;
};

Shader *Fur2Phase::createShader(Renderer *renderer) const {
	return new Fur2Shader(renderer, colorTexture.get());
}

MTS_IMPLEMENT_CLASS(Fur2Shader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Fur2Phase, false, SpectrumPhaseFunction)
MTS_EXPORT_PLUGIN(Fur2Phase, "Our fur2 BSDF")
MTS_NAMESPACE_END
