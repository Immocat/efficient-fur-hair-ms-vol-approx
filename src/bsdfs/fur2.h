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

#include "fur2_basics.h"
#include "mitsuba/render/bsdf.h"
#include "mitsuba/render/scene.h"
#include "mitsuba/render/texture.h"

MTS_NAMESPACE_BEGIN

#ifdef _WIN32
#ifdef FUR2_EXPORTS
#define FUR2_API __declspec(dllexport)
#else
#define FUR2_API __declspec(dllimport)
#endif
#else
#define FUR2_API __attribute__((visibility("default")))
#endif

class FUR2_API Fur2 : public BSDF {
  public:
	Fur2(const Properties &props);

	Fur2(Stream *stream, InstanceManager *manager);

	void configure();

	Spectrum getDiffuseReflectance(const Intersection &its) const;

	Vector3 normalizeSafe(Vector3 v) const;

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const;

	Float normrnd(const BSDFSamplingRecord &bRec, Float mu, Float sigma) const;

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const;

	// Choose a lobe and sample it.
	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const;

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const;

	Spectrum sampleWithWorldWiAndOrientation(const Vector3 &worldWi, const Vector3 &orientation, Float &pdf, Vector3f &worldWo, Sampler *sampler) const;

	Spectrum sampleWithLocalWi(const Vector3 &wi, Vector3 &wo, Float &pdf, const Point2 &sample, Float &eta, BSDF::EBSDFType &sampledType, Sampler *sampler) const;

	void addChild(const std::string &name, ConfigurableObject *child);

	void serialize(Stream *stream, InstanceManager *manager) const;

	Float getRoughness(const Intersection &its, int component) const;

	std::string toString() const;

	Shader *createShader(Renderer *renderer) const;

	MatrixRGB PrecomputeHairMatrixRGB(int bands, int nSamples) const;

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
};

MTS_NAMESPACE_END