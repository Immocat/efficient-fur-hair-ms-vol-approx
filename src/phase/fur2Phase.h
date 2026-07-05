#pragma once

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

#include "mitsuba/render/bsdf.h"
#include "mitsuba/render/scene.h"
#include "mitsuba/render/texture.h"
#include "src/bsdfs/fur2_basics.h"

MTS_NAMESPACE_BEGIN

class Fur2Phase : public SpectrumPhaseFunction {
  public:
	Fur2Phase(const Properties &props);

	Fur2Phase(Stream *stream, InstanceManager *manager);

	void configure();

	Spectrum getDiffuseReflectance(const Intersection &its) const;

	Vector3 normalizeSafe(Vector3 v) const;

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const;

	Float normrnd(const BSDFSamplingRecord &bRec, Float mu, Float sigma) const;

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const;

	// Choose a lobe and sample it.
	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const;

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const;

	/// Phase function interfaces
	Float eval(const PhaseFunctionSamplingRecord &pRec) const;

	Float sample(PhaseFunctionSamplingRecord &pRec,
				 Sampler *sampler) const;

	Float sample(PhaseFunctionSamplingRecord &pRec,
				 Float &pdf, Sampler *sampler) const;

	Vector3 evalRGB(const PhaseFunctionSamplingRecord &pRec, const Float sample1D) const override;

	Vector3 lookupTotalScatteredEnergy(Float thetaI) const;

	bool needsDirectionallyVaryingCoefficients() const override;
	Float sigmaDir(Float cosTheta) const override;
	Float sigmaDirMax() const override;
	Float getMeanCosine() const override;
	Spectrum albedoDir(const Vector &d, const Vector &orientation) const override;
	bool needsDirectionallyVaryingAlbedo() const override;

	void addChild(const std::string &name, ConfigurableObject *child);

	void serialize(Stream *stream, InstanceManager *manager) const;

	Float getRoughness(const Intersection &its, int component) const;

	std::string toString() const;

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
  private:
	std::string renderMode;
	Float colorScale, colorGamma;
	ref<Texture> colorTexture;
	bool useColorTexture;
	Float medulla_ratio;
	Float eta;
	Vector3 absorption_outer;
	Float beta_m, beta_n;
	Float scattering_inner, absorption_inner;
	Float alpha;
	Float g;
	Float layers;
	Float absorptionFactor;
	int furMaskCode;
	PerspectiveCamera *perspectiveCamera;

	PFPrecalcSample integrateBSDF(const Vector3f &wi, Sampler *sampler) const;
	void generateIntegratedHairBSDFTable();
	PFPrecalcSample sampleWithLocalWiForIntegration(const Vector3 &wi, Vector3 &wo, Float &pdf, Sampler *sampler) const;

	/// Phase
	int m_theta_i_steps = 480.f;
	int m_intergration_sample_cnt = 65536;

	/// Storage for the precomputed table
	std::vector<Float> m_normalizationFactors;
	std::shared_ptr<std::vector<Vector3>> m_totalScatteredEnergy;
};

MTS_NAMESPACE_END