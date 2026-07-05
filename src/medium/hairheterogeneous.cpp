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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/statistics.h>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

MTS_NAMESPACE_BEGIN

#define HETVOL_EARLY_EXIT 1

class HairHeterogeneousMedium : public Medium {
public:
	enum EIntegrationMethod {
		ESimpsonQuadrature = 0,
		EWoodcockTracking
	};

	HairHeterogeneousMedium(const Properties &props)
		: Medium(props) {
		m_stepSize = props.getFloat("stepSize", 0);
		m_scale = props.getFloat("scale", 1);
		m_useTensorDirectionalDensity = props.getBoolean("useTensorDirectionalDensity", false);

		if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
			Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
						"homogeneous media. Please use nested volume instances.");
		if (props.hasProperty("densityMultiplier"))
			Log(EError, "The 'densityMultiplier' parameter is deprecated and now called 'scale'.");

		std::string method = boost::to_lower_copy(props.getString("method", "woodcock"));
		if (method == "woodcock")
			m_method = EWoodcockTracking;
		else if (method == "simpson")
			m_method = ESimpsonQuadrature;
		else
			Log(EError, "Unsupported integration method \"%s\"!", method.c_str());
	}

	HairHeterogeneousMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager) {
		m_method = (EIntegrationMethod)stream->readInt();
		m_scale = stream->readFloat();
		m_density = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_orientation = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_variance = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_useTensorDirectionalDensity = stream->readBool();
		m_tensorDiag = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_tensorOffDiag = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		m_stepSize = stream->readFloat();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
		stream->writeInt(m_method);
		stream->writeFloat(m_scale);
		manager->serialize(stream, m_density.get());
		manager->serialize(stream, m_albedo.get());
		manager->serialize(stream, m_orientation.get());
		manager->serialize(stream, m_variance.get());
		stream->writeBool(m_useTensorDirectionalDensity);
		manager->serialize(stream, m_tensorDiag.get());
		manager->serialize(stream, m_tensorOffDiag.get());
		stream->writeFloat(m_stepSize);
	}

	void configure() {
		Medium::configure();
		if (m_density.get() == NULL)
			Log(EError, "No density specified!");
		if (m_albedo.get() == NULL)
			Log(EError, "No albedo specified!");

		m_densityAABB = m_density->getAABB();
		m_anisotropicMedium = m_phaseFunction->needsDirectionallyVaryingCoefficients();
		m_anisotropicAlbedo = m_phaseFunction->needsDirectionallyVaryingAlbedo();

		m_maxDensity = m_scale * m_density->getMaximumFloatValue();
		if (m_anisotropicMedium && !m_useTensorDirectionalDensity)
			m_maxDensity *= m_phaseFunction->sigmaDirMax();
		m_invMaxDensity = 1.0f / m_maxDensity;

		if (m_stepSize == 0) {
			m_stepSize = std::min(m_density->getStepSize(), m_albedo->getStepSize());
			if (m_orientation != NULL)
				m_stepSize = std::min(m_stepSize, m_orientation->getStepSize());
			if (m_tensorDiag != NULL)
				m_stepSize = std::min(m_stepSize, m_tensorDiag->getStepSize());
			if (m_tensorOffDiag != NULL)
				m_stepSize = std::min(m_stepSize, m_tensorOffDiag->getStepSize());

			if (m_stepSize == std::numeric_limits<Float>::infinity())
				Log(EError, "Unable to infer a suitable step size. Please set 'stepSize'.");
		}

		if (m_useTensorDirectionalDensity) {
			if (m_tensorDiag.get() == NULL || m_tensorOffDiag.get() == NULL)
				Log(EError, "Tensor directional density requires both 'tensorDiag' and 'tensorOffDiag' volumes.");
		} else {
			if (m_anisotropicMedium && m_orientation.get() == NULL)
				Log(EError, "Cannot use anisotropic phase function without orientation volume.");
			if (m_anisotropicAlbedo && m_orientation.get() == NULL)
				Log(EError, "Cannot use anisotropic albedo without orientation volume.");
		}
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			VolumeDataSource *volume = static_cast<VolumeDataSource *>(child);
			if (name == "albedo") {
				Assert(volume->supportsSpectrumLookups());
				m_albedo = volume;
			} else if (name == "density") {
				Assert(volume->supportsFloatLookups());
				m_density = volume;
			} else if (name == "orientation") {
				Assert(volume->supportsVectorLookups());
				m_orientation = volume;
			} else if (name == "variance") {
				Assert(volume->supportsFloatLookups());
				m_variance = volume;
			} else if (name == "tensorDiag") {
				Assert(volume->supportsVectorLookups());
				m_tensorDiag = volume;
			} else if (name == "tensorOffDiag") {
				Assert(volume->supportsVectorLookups());
				m_tensorOffDiag = volume;
			} else {
				Medium::addChild(name, child);
			}
		} else {
			Medium::addChild(name, child);
		}
	}

	Float integrateDensity(const Ray &ray) const {
		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return 0.0f;

		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt - mint, maxComp = 0;
		Point p = ray(mint), pLast = ray(maxt);

		for (int i = 0; i < 3; ++i)
			maxComp = std::max(std::max(maxComp, std::abs(p[i])), std::abs(pLast[i]));
		if (length < 1e-6f * maxComp)
			return 0.0f;

		uint32_t nSteps = (uint32_t)std::ceil(length / m_stepSize);
		nSteps += nSteps % 2;
		const Float stepSize = length / nSteps;
		const Vector increment = ray.d * stepSize;

		Float integratedDensity = lookupDensity(p, ray.d) + lookupDensity(pLast, ray.d);

#if defined(HETVOL_EARLY_EXIT)
		const Float stopAfterDensity = -math::fastlog(Epsilon);
		const Float stopValue = stopAfterDensity * 3.0f / (stepSize * m_scale);
#endif

		p += increment;
		Float m = 4;
		for (uint32_t i = 1; i < nSteps; ++i) {
			integratedDensity += m * lookupDensity(p, ray.d);
			m = 6 - m;

#if defined(HETVOL_EARLY_EXIT)
			if (integratedDensity > stopValue)
				return std::numeric_limits<Float>::infinity();
#endif

			Point next = p + increment;
			if (p == next) {
				Log(EWarn, "integrateDensity(): unable to make forward progress.");
				break;
			}
			p = next;
		}

		return integratedDensity * m_scale * stepSize * (1.0f / 3.0f);
	}

	bool invertDensityIntegral(const Ray &ray, Float desiredDensity,
							   Float &integratedDensity, Float &t, Float &densityAtMinT,
							   Float &densityAtT) const {
		integratedDensity = densityAtMinT = densityAtT = 0.0f;

		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return false;
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);
		Float length = maxt - mint, maxComp = 0;
		Point p = ray(mint), pLast = ray(maxt);

		for (int i = 0; i < 3; ++i)
			maxComp = std::max(std::max(maxComp, std::abs(p[i])), std::abs(pLast[i]));
		if (length < 1e-6f * maxComp)
			return 0.0f;

		uint32_t nSteps = (uint32_t)std::ceil(length / (2 * m_stepSize));
		Float stepSize = length / nSteps, multiplier = (1.0f / 6.0f) * stepSize * m_scale;
		Vector fullStep = ray.d * stepSize, halfStep = fullStep * .5f;

		Float node1 = lookupDensity(p, ray.d);
		if (ray.mint == mint)
			densityAtMinT = node1 * m_scale;
		else
			densityAtMinT = 0.0f;

		for (uint32_t i = 0; i < nSteps; ++i) {
			Float node2 = lookupDensity(p + halfStep, ray.d),
				  node3 = lookupDensity(p + fullStep, ray.d),
				  newDensity = integratedDensity + multiplier * (node1 + node2 * 4 + node3);

			if (newDensity >= desiredDensity) {
				Float a = 0, b = stepSize, x = a, fx = integratedDensity - desiredDensity,
					  stepSizeSqr = stepSize * stepSize, temp = m_scale / stepSizeSqr;
				int it = 1;

				while (true) {
					Float dfx = temp * (node1 * stepSizeSqr -
										(3 * node1 - 4 * node2 + node3) * stepSize * x +
										2 * (node1 - 2 * node2 + node3) * x * x);
					x -= fx / dfx;

					if (EXPECT_NOT_TAKEN(x <= a || x >= b || dfx == 0))
						x = 0.5f * (b + a);

					Float intval = integratedDensity + temp * (1.0f / 6.0f) * (x *
																				 (6 * node1 * stepSizeSqr -
																				  3 * (3 * node1 - 4 * node2 + node3) * stepSize * x +
																				  4 * (node1 - 2 * node2 + node3) * x * x));
					fx = intval - desiredDensity;

					if (std::abs(fx) < 1e-6f) {
						t = mint + stepSize * i + x;
						integratedDensity = intval;
						densityAtT = temp * (node1 * stepSizeSqr -
											 (3 * node1 - 4 * node2 + node3) * stepSize * x +
											 2 * (node1 - 2 * node2 + node3) * x * x);
						return true;
					} else if (++it > 30) {
						Log(EWarn, "invertDensityIntegral(): Newton-Bisection failed.");
						return false;
					}

					if (fx > 0)
						b = x;
					else
						a = x;
				}
			}

			Point next = p + fullStep;
			if (p == next) {
				Log(EWarn, "invertDensityIntegral(): unable to make forward progress.");
				break;
			}
			integratedDensity = newDensity;
			node1 = node3;
			p = next;
		}
		return false;
	}

	Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const {
		if (m_method == ESimpsonQuadrature || sampler == NULL) {
			return Spectrum(math::fastexp(-integrateDensity(ray)));
		}

		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return Spectrum(1.0f);
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);

		int nSamples = 2;
		Float result = 0;
		for (int i = 0; i < nSamples; ++i) {
			Float tDummy, densityDummy;
			Float acceptDummy = 0.0f, maxRatioDummy = 0.0f;
			bool collided = m_density->supportsFloatLocalMajorant()
								? findCollisionWoodcockLocal(ray, sampler, tDummy, densityDummy, acceptDummy, maxRatioDummy, mint, maxt)
								: findCollisionWoodcockGlobal(ray, sampler, tDummy, densityDummy, acceptDummy, maxRatioDummy, mint, maxt);
			if (!collided)
				result += 1;
		}
		return Spectrum(result / nSamples);
	}

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec, Sampler *sampler) const {
		Float integratedDensity, densityAtMinT, densityAtT;
		bool success = false;
		mRec.variance = 0.0f;
		mRec.profileAcceptRatio = 0.0f;
		mRec.profileMaxRatio = 0.0f;

		if (m_method == ESimpsonQuadrature) {
			Float desiredDensity = -math::fastlog(1 - sampler->next1D());
			if (invertDensityIntegral(ray, desiredDensity, integratedDensity,
									  mRec.t, densityAtMinT, densityAtT)) {
				mRec.p = ray(mRec.t);
				success = true;
				Spectrum albedo = lookupAlbedo(mRec.p, ray.d);
				mRec.sigmaS = albedo * densityAtT;
				mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
				mRec.orientation = getOrientation(mRec.p);
				mRec.variance = m_variance != nullptr ? m_variance->lookupFloat(mRec.p) : Float(0.f);
			}

			Float expVal = math::fastexp(-integratedDensity);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * densityAtT;
			mRec.pdfSuccessRev = expVal * densityAtMinT;
			mRec.transmittance = Spectrum(expVal);
			mRec.time = ray.time;
			mRec.profileAcceptRatio = 0.0f;
			mRec.profileMaxRatio = 0.0f;
		} else {
			mRec.pdfFailure = 1.0f;
			mRec.pdfSuccess = 1.0f;
			mRec.pdfSuccessRev = 1.0f;
			mRec.transmittance = Spectrum(1.0f);
			mRec.time = ray.time;
			mRec.profileAcceptRatio = 0.0f;
			mRec.profileMaxRatio = 0.0f;

			Float mint, maxt;
			if (!m_densityAABB.rayIntersect(ray, mint, maxt))
				return false;
			mint = std::max(mint, ray.mint);
			maxt = std::min(maxt, ray.maxt);

			Float hitT = 0, hitDensity = 0;
			Float hitAcceptRatio = 0, maxRatioAlongMarch = 0;
			success = m_density->supportsFloatLocalMajorant()
						  ? findCollisionWoodcockLocal(ray, sampler, hitT, hitDensity, hitAcceptRatio, maxRatioAlongMarch, mint, maxt)
						  : findCollisionWoodcockGlobal(ray, sampler, hitT, hitDensity, hitAcceptRatio, maxRatioAlongMarch, mint, maxt);
			if (success) {
				mRec.t = hitT;
				mRec.p = ray(hitT);
				Spectrum albedo = lookupAlbedo(mRec.p, ray.d);
				mRec.sigmaS = albedo * hitDensity;
				mRec.sigmaA = Spectrum(hitDensity) - mRec.sigmaS;
				mRec.transmittance = Spectrum(hitDensity != 0.0f ? 1.0f / hitDensity : 0);
				if (!std::isfinite(mRec.transmittance[0]))
					mRec.transmittance = Spectrum(0.0f);
				mRec.orientation = getOrientation(mRec.p);
				mRec.variance = m_variance != nullptr ? m_variance->lookupFloat(mRec.p) : Float(0.f);
				mRec.medium = this;
			}
			mRec.profileAcceptRatio = hitAcceptRatio;
			mRec.profileMaxRatio = maxRatioAlongMarch;
		}
		mRec.medium = this;
		return success && mRec.pdfSuccess > 0;
	}

	void eval(const Ray &ray, MediumSamplingRecord &mRec) const {
		if (m_method == ESimpsonQuadrature) {
			Float expVal = math::fastexp(-integrateDensity(ray));
			Float mintDensity = lookupDensity(ray(ray.mint), ray.d) * m_scale;
			Float maxtDensity = 0.0f;
			Spectrum maxtAlbedo(0.0f);
			if (ray.maxt < std::numeric_limits<Float>::infinity()) {
				Point p = ray(ray.maxt);
				maxtDensity = lookupDensity(p, ray.d) * m_scale;
				maxtAlbedo = lookupAlbedo(p, ray.d);
			}
			mRec.transmittance = Spectrum(expVal);
			mRec.pdfFailure = expVal;
			mRec.pdfSuccess = expVal * maxtDensity;
			mRec.pdfSuccessRev = expVal * mintDensity;
			mRec.sigmaS = maxtAlbedo * maxtDensity;
			mRec.sigmaA = Spectrum(maxtDensity) - mRec.sigmaS;
			mRec.time = ray.time;
			mRec.profileAcceptRatio = 0.0f;
			mRec.profileMaxRatio = 0.0f;
			mRec.medium = this;
		} else {
			Log(EError, "eval(): unsupported integration method!");
		}
	}

	bool isHomogeneous() const {
		return false;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HairHeterogeneousMedium[" << endl
			<< "  density = " << indent(m_density.toString()) << "," << endl
			<< "  albedo = " << indent(m_albedo.toString()) << "," << endl
			<< "  orientation = " << indent(m_orientation.toString()) << "," << endl
			<< "  tensorDiag = " << indent(m_tensorDiag.toString()) << "," << endl
			<< "  tensorOffDiag = " << indent(m_tensorOffDiag.toString()) << "," << endl
			<< "  useTensorDirectionalDensity = " << m_useTensorDirectionalDensity << "," << endl
			<< "  stepSize = " << m_stepSize << "," << endl
			<< "  scale = " << m_scale << endl
			<< "]";
		return oss.str();
	}

	void evalCoefficientsWithDirection(const Point &p, const Vector3 &d,
									   Spectrum &sigmaS, Spectrum &sigmaT, Vector3 &orientation) const {
		Float density = lookupDensity(p, d) * m_scale;
		orientation = getOrientation(p);

		if (density < Epsilon) {
			sigmaS = Spectrum(0.0f);
			sigmaT = Spectrum(0.0f);
			return;
		}

		Spectrum albedo = lookupAlbedo(p, d);
		sigmaS = albedo * density;
		sigmaT = Spectrum(density);
	}

	void evalCoefficients(const Point &p, Spectrum &sigmaS, Spectrum &sigmaT) const {
		Float density = m_density->lookupFloat(p) * m_scale;
		if (density < Epsilon) {
			sigmaS = Spectrum(0.0f);
			sigmaT = Spectrum(0.0f);
			return;
		}

		Spectrum albedo = m_albedo->lookupSpectrum(p);
		sigmaS = albedo * density;
		sigmaT = Spectrum(density);
	}

	void evalDensityWithDirection(const Point &p, const Vector3 &d, Float &sigmaT) const {
		sigmaT = lookupDensity(p, d) * m_scale;
	}

	Vector3 getOrientation(const Point &p) const {
		if (m_useTensorDirectionalDensity)
			return principalDirectionFromTensor(p);
		if (m_orientation == NULL)
			return Vector3(0.0f);
		return m_orientation->lookupVector(p);
	}

	Float getOrientationVariance(const Point &p) const {
		if (m_variance == nullptr) {
			Log(EError, "m_variance is null!");
			return Float(0.f);
		}
		return m_variance->lookupFloat(p);
	}

	MTS_DECLARE_CLASS()
protected:
	inline bool lookupTensorAt(const Point &p, Vector3 &diag, Vector3 &offDiag) const {
		if (m_tensorDiag == NULL || m_tensorOffDiag == NULL)
			return false;
		diag = m_tensorDiag->lookupVector(p);
		offDiag = m_tensorOffDiag->lookupVector(p);
		return true;
	}

	inline Vector3 mulSymmetricTensor(const Vector3 &diag, const Vector3 &offDiag, const Vector3 &v) const {
		return Vector3(
			diag.x * v.x + offDiag.x * v.y + offDiag.y * v.z,
			offDiag.x * v.x + diag.y * v.y + offDiag.z * v.z,
			offDiag.y * v.x + offDiag.z * v.y + diag.z * v.z);
	}

	inline Vector3 principalDirectionFromTensor(const Point &p) const {
		Vector3 diag, offDiag;
		if (!lookupTensorAt(p, diag, offDiag))
			return Vector3(0.0f);

		Vector3 v(1.0f, 0.0f, 0.0f);
		for (int i = 0; i < 8; ++i) {
			Vector3 Tv = mulSymmetricTensor(diag, offDiag, v);
			Float len = Tv.length();
			v = (len > 1e-6f) ? (Tv / len) : Vector3(1.0f, 0.0f, 0.0f);
		}
		return v.lengthSquared() > 1e-12f ? normalize(v) : Vector3(0.0f);
	}

	inline Float tensorDirectionalFactor(const Point &p, const Vector &d) const {
		Vector3 diag, offDiag;
		if (!lookupTensorAt(p, diag, offDiag))
			return 0.0f;

		Vector w = normalize(d);
		Float pParallel =
			diag.x * w.x * w.x +
			diag.y * w.y * w.y +
			diag.z * w.z * w.z +
			2.0f * offDiag.x * w.x * w.y +
			2.0f * offDiag.y * w.x * w.z +
			2.0f * offDiag.z * w.y * w.z;
		pParallel = math::clamp(pParallel, (Float)0.0f, (Float)1.0f);
		return std::sqrt(std::max((Float)0.0f, (Float)1.0f - pParallel));
	}

	inline Float lookupDensity(const Point &p, const Vector &d) const {
		Float density = m_density->lookupFloat(p);
		if (m_useTensorDirectionalDensity) {
			if (density == 0.0f)
				return 0.0f;
			return density * tensorDirectionalFactor(p, d);
		}

		if (m_anisotropicMedium && density != 0) {
			Vector orientation = m_orientation->lookupVector(p);
			if (!orientation.isZero())
				density *= m_phaseFunction->sigmaDir(dot(d, orientation));
			else
				return 0;
		}
		return density;
	}

	inline Spectrum lookupAlbedo(const Point &p, const Vector &d) const {
		Spectrum albedo = m_albedo->lookupSpectrum(p);
		if (m_anisotropicAlbedo && !albedo.isZero()) {
			Vector orientation = getOrientation(p);
			if (!orientation.isZero())
				albedo *= m_phaseFunction->albedoDir(d, orientation);
			else
				return Spectrum(0.0f);
		}
		return albedo;
	}

	bool findCollisionWoodcockGlobal(const Ray &ray, Sampler *sampler,
									 Float &hitT, Float &hitDensity,
									 Float &hitAcceptRatio, Float &maxRatioAlongMarch,
									 Float mint, Float maxt) const {
		Float t = mint;
		hitAcceptRatio = 0.0f;
		maxRatioAlongMarch = 0.0f;
		while (true) {
			t -= math::fastlog(1 - sampler->next1D()) * m_invMaxDensity;
			if (t >= maxt)
				break;

			Point p = ray(t);
			Float densityAtT = lookupDensity(p, ray.d) * m_scale;
			Float ratio = densityAtT * m_invMaxDensity;
			maxRatioAlongMarch = std::max(maxRatioAlongMarch, ratio);
			if (ratio > sampler->next1D()) {
				hitT = t;
				hitDensity = densityAtT;
				hitAcceptRatio = ratio;
				return true;
			}
		}
		return false;
	}

	bool findCollisionWoodcockLocal(const Ray &ray, Sampler *sampler,
									Float &hitT, Float &hitDensity,
									Float &hitAcceptRatio, Float &maxRatioAlongMarch,
									Float mint, Float maxt) const {
		Float t = mint;
		const Float epsAdvance = std::max((Float)1e-7f, (Float)1e-6f * (maxt - mint));
		const Float anisotropicBound = (m_anisotropicMedium && !m_useTensorDirectionalDensity) ? m_phaseFunction->sigmaDirMax() : 1.0f;
		hitAcceptRatio = 0.0f;
		maxRatioAlongMarch = 0.0f;

		while (t < maxt) {
			const Float minAdvance = std::max((Float)1e-9f,
											  std::abs(t) * std::numeric_limits<Float>::epsilon() * (Float)8.0f);
			Point p = ray(t);
			Float localMajorant, segLen;
			if (!m_density->lookupFloatLocalMajorantSegment(p, ray.d, localMajorant, segLen) || segLen <= 0) {
				t += epsAdvance;
				continue;
			}

			Float majorant = localMajorant * m_scale * anisotropicBound;
			Float tSegEnd = std::min(maxt, t + segLen);
			if (tSegEnd <= t + minAdvance) {
				t += epsAdvance;
				continue;
			}
			if (majorant <= 0) {
				t = tSegEnd + epsAdvance;
				continue;
			}

			while (t < tSegEnd) {
				Float u = sampler->next1D();
				u = std::max((Float)0.0f, std::min((Float)(1.0f - 1e-7f), u));
				Float step = -math::fastlog(1 - u) / majorant;
				if (step <= 0)
					step = minAdvance;
				Float tCandidate = t + step;
				if (tCandidate <= t + minAdvance)
					tCandidate = t + minAdvance;
				if (tCandidate >= tSegEnd)
					break;

				Float densityAtT = lookupDensity(ray(tCandidate), ray.d) * m_scale;
				Float ratio = densityAtT / majorant;
				maxRatioAlongMarch = std::max(maxRatioAlongMarch, ratio);
				if (ratio > sampler->next1D()) {
					hitT = tCandidate;
					hitDensity = densityAtT;
					hitAcceptRatio = ratio;
					return true;
				}
				t = tCandidate;
			}

			t = (tSegEnd > t + minAdvance) ? tSegEnd : t + epsAdvance;
		}
		return false;
	}

protected:
	EIntegrationMethod m_method;
	ref<VolumeDataSource> m_density;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientation;
	ref<VolumeDataSource> m_variance;
	ref<VolumeDataSource> m_tensorDiag;
	ref<VolumeDataSource> m_tensorOffDiag;
	Float m_scale;
	bool m_anisotropicMedium;
	bool m_anisotropicAlbedo;
	bool m_useTensorDirectionalDensity;
	Float m_stepSize;
	AABB m_densityAABB;
	Float m_maxDensity;
	Float m_invMaxDensity;
};

MTS_IMPLEMENT_CLASS_S(HairHeterogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HairHeterogeneousMedium, "Hair heterogeneous medium with local-majorant Woodcock")
MTS_NAMESPACE_END
