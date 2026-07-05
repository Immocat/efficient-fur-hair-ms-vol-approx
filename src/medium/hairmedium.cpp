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

#include "hairmedium.h"
#include <boost/algorithm/string.hpp>
#include <mitsuba/core/statistics.h>
#include <mitsuba/render/scene.h>
// #pragma optimize("", off)
MTS_NAMESPACE_BEGIN

/**
 * \brief When the following line is uncommented, the medium implementation
 * stops integrating density when it is determined that the segment has a
 * throughput of less than 'Epsilon' (see \c mitsuba/core/constants.h)
 */
#define HETVOL_EARLY_EXIT 1

/// Generate a few statistics related to the implementation?
// #define HETVOL_STATISTICS 1

#if defined(HETVOL_STATISTICS)
static StatsCounter avgNewtonIterations("Heterogeneous volume",
										"Avg. # of Newton-Bisection iterations", EAverage);
static StatsCounter avgRayMarchingStepsTransmittance("Heterogeneous volume",
													 "Avg. # of ray marching steps (transmittance)", EAverage);
static StatsCounter avgRayMarchingStepsSampling("Heterogeneous volume",
												"Avg. # of ray marching steps (sampling)", EAverage);
static StatsCounter earlyExits("Heterogeneous volume",
							   "Number of early exits", EPercentage);
#endif

HairMedium::HairMedium(const Properties &props)
	: Medium(props) {
	m_stepSize = props.getFloat("stepSize", 0);
	m_scale = props.getFloat("scale", 1);
	if (props.hasProperty("sigmaS") || props.hasProperty("sigmaA"))
		Log(EError, "The 'sigmaS' and 'sigmaA' properties are only supported by "
					"homogeneous media. Please use nested volume instances to supply "
					"these parameters");

	if (props.hasProperty("densityMultiplier"))
		Log(EError, "The 'densityMultiplier' parameter has been deprecated and is now called 'scale'.");

	std::string method = boost::to_lower_copy(props.getString("method", "woodcock"));
	if (method == "woodcock")
		m_method = EWoodcockTracking;
	else if (method == "simpson")
		m_method = ESimpsonQuadrature;
	else
		Log(EError, "Unsupported integration method \"%s\"!", method.c_str());

	m_LPVRes = props.getInteger("LPVResolution", 8);
	m_LPV = VoxelGrid<SHVectorRGB>(Vector3i(m_LPVRes));
	m_biasedStepping = props.getBoolean("biasedStepping", false);
}

/* Unserialize from a binary data stream */
HairMedium::HairMedium(Stream *stream, InstanceManager *manager)
	: Medium(stream, manager) {
	m_method = (EIntegrationMethod)stream->readInt();
	m_scale = stream->readFloat();
	m_density = static_cast<VolumeDataSource *>(manager->getInstance(stream));
	m_albedo = static_cast<VolumeDataSource *>(manager->getInstance(stream));
	m_orientation = static_cast<VolumeDataSource *>(manager->getInstance(stream));
	m_stepSize = stream->readFloat();
	configure();
}

HairMedium::~HairMedium() {
}

/* Serialize the volume to a binary data stream */
void HairMedium::serialize(Stream *stream, InstanceManager *manager) const {
	Medium::serialize(stream, manager);
	stream->writeInt(m_method);
	stream->writeFloat(m_scale);
	manager->serialize(stream, m_density.get());
	manager->serialize(stream, m_albedo.get());
	manager->serialize(stream, m_orientation.get());
	stream->writeFloat(m_stepSize);
}

void HairMedium::configure() {
	Medium::configure();
	if (m_density.get() == NULL)
		Log(EError, "No density specified!");
	if (m_albedo.get() == NULL)
		Log(EError, "No albedo specified!");
	m_densityAABB = m_density->getAABB();
	m_anisotropicMedium =
		m_phaseFunction->needsDirectionallyVaryingCoefficients();
	m_anisotropicAlbedo =
		m_phaseFunction->needsDirectionallyVaryingAlbedo();

	/* Assumes that the density medium does not
	   contain values greater than one! */
	m_maxDensity = m_scale * m_density->getMaximumFloatValue();
	if (m_anisotropicMedium)
		m_maxDensity *= m_phaseFunction->sigmaDirMax();
	m_invMaxDensity = 1.0f / m_maxDensity;
	if (m_biasedStepping) {
		m_invMaxDensity = 1.f;
	}
	if (m_stepSize == 0) {
		m_stepSize = std::min(
			m_density->getStepSize(), m_albedo->getStepSize());
		if (m_orientation != NULL)
			m_stepSize = std::min(m_stepSize,
								  m_orientation->getStepSize());

		if (m_stepSize == std::numeric_limits<Float>::infinity())
			Log(EError, "Unable to infer a suitable step size for deterministic "
						"integration, please specify one manually using the 'stepSize' "
						"parameter.");
	}

	if (m_anisotropicMedium && m_orientation.get() == NULL)
		Log(EError, "Cannot use anisotropic phase function: "
					"did not specify a particle orientation field!");

	if (m_anisotropicAlbedo && m_orientation.get() == NULL)
		Log(EError, "Cannot use anisotropic albedo: "
					"did not specify a particle orientation field!");
}

void HairMedium::addChild(const std::string &name, ConfigurableObject *child) {
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
		} else {
			Medium::addChild(name, child);
		}
	} else {
		Medium::addChild(name, child);
	}
}

inline Float HairMedium::lookupDensity(const Point &p, const Vector &d) const {
	Float density = m_density->lookupFloat(p);
	if (m_anisotropicMedium && density != 0) {
		Vector orientation = m_orientation->lookupVector(p);
		if (!orientation.isZero())
			density *= m_phaseFunction->sigmaDir(dot(d, orientation));
		else
			return 0;
	}
	return density;
}

inline Spectrum HairMedium::lookupAlbedo(const Point &p, const Vector &d) const {
	Spectrum albedo = m_albedo->lookupSpectrum(p);
	if (m_anisotropicAlbedo && !albedo.isZero()) {
		Vector orientation = m_orientation->lookupVector(p);
		if (!orientation.isZero())
			albedo *= m_phaseFunction->albedoDir(d, orientation);
		else
			return Spectrum(0.0f);
	}
	return albedo;
}

Float HairMedium::integrateDensity(const Ray &ray) const {
	/* Determine the ray segment, along which the
	   density integration should take place */
	Float mint, maxt;
	if (!m_densityAABB.rayIntersect(ray, mint, maxt))
		return 0.0f;

	mint = std::max(mint, ray.mint);
	maxt = std::min(maxt, ray.maxt);
	Float length = maxt - mint, maxComp = 0;

	Point p = ray(mint), pLast = ray(maxt);

	/* Ignore degenerate path segments */
	for (int i = 0; i < 3; ++i)
		maxComp = std::max(std::max(maxComp,
									std::abs(p[i])),
						   std::abs(pLast[i]));
	if (length < 1e-6f * maxComp)
		return 0.0f;

	/* Compute a suitable step size */
	uint32_t nSteps = (uint32_t)std::ceil(length / m_stepSize);
	nSteps += nSteps % 2;
	const Float stepSize = length / nSteps;
	const Vector increment = ray.d * stepSize;

#if defined(HETVOL_STATISTICS)
	avgRayMarchingStepsTransmittance.incrementBase();
	earlyExits.incrementBase();
#endif

	/* Perform lookups at the first and last node */
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

#if defined(HETVOL_STATISTICS)
		++avgRayMarchingStepsTransmittance;
#endif

#if defined(HETVOL_EARLY_EXIT)
		if (integratedDensity > stopValue) {
// Reached the threshold -- stop early
#if defined(HETVOL_STATISTICS)
			++earlyExits;
#endif
			return std::numeric_limits<Float>::infinity();
		}
#endif

		Point next = p + increment;
		if (p == next) {
			Log(EWarn, "integrateDensity(): unable to make forward progress -- "
					   "round-off error issues? The step size was %e, mint=%f, "
					   "maxt=%f, nSteps=%i, ray=%s",
				stepSize, mint, maxt, nSteps,
				ray.toString().c_str());
			break;
		}
		p = next;
	}

	return integratedDensity * m_scale * stepSize * (1.0f / 3.0f);
}

bool HairMedium::invertDensityIntegral(const Ray &ray, Float desiredDensity,
									   Float &integratedDensity, Float &t, Float &densityAtMinT,
									   Float &densityAtT) const {
	integratedDensity = densityAtMinT = densityAtT = 0.0f;

	/* Determine the ray segment, along which the
	   density integration should take place */
	Float mint, maxt;
	if (!m_densityAABB.rayIntersect(ray, mint, maxt))
		return false;
	mint = std::max(mint, ray.mint);
	maxt = std::min(maxt, ray.maxt);
	Float length = maxt - mint, maxComp = 0;
	Point p = ray(mint), pLast = ray(maxt);

	/* Ignore degenerate path segments */
	for (int i = 0; i < 3; ++i)
		maxComp = std::max(std::max(maxComp,
									std::abs(p[i])),
						   std::abs(pLast[i]));
	if (length < 1e-6f * maxComp)
		return 0.0f;

	/* Compute a suitable step size (this routine samples the integrand
	   between steps, hence the factor of 2) */
	uint32_t nSteps = (uint32_t)std::ceil(length / (2 * m_stepSize));
	Float stepSize = length / nSteps,
		  multiplier = (1.0f / 6.0f) * stepSize * m_scale;
	Vector fullStep = ray.d * stepSize,
		   halfStep = fullStep * .5f;

	Float node1 = lookupDensity(p, ray.d);

	if (ray.mint == mint)
		densityAtMinT = node1 * m_scale;
	else
		densityAtMinT = 0.0f;

#if defined(HETVOL_STATISTICS)
	avgRayMarchingStepsSampling.incrementBase();
#endif

	for (uint32_t i = 0; i < nSteps; ++i) {
		Float node2 = lookupDensity(p + halfStep, ray.d),
			  node3 = lookupDensity(p + fullStep, ray.d),
			  newDensity = integratedDensity + multiplier *
												   (node1 + node2 * 4 + node3);
#if defined(HETVOL_STATISTICS)
		++avgRayMarchingStepsSampling;
#endif
		if (newDensity >= desiredDensity) {
			/* The integrated density of the last segment exceeds the desired
			   amount -- now use the Simpson quadrature expression and
			   Newton-Bisection to find the precise location of the scattering
			   event. Note that no further density queries are performed after
			   this point; instead, the density are modeled based on a
			   quadratic polynomial that is fit to the last three lookups */

			Float a = 0, b = stepSize, x = a,
				  fx = integratedDensity - desiredDensity,
				  stepSizeSqr = stepSize * stepSize,
				  temp = m_scale / stepSizeSqr;
			int it = 1;

#if defined(HETVOL_STATISTICS)
			avgNewtonIterations.incrementBase();
#endif
			while (true) {
#if defined(HETVOL_STATISTICS)
				++avgNewtonIterations;
#endif
				/* Lagrange polynomial from the Simpson quadrature */
				Float dfx = temp * (node1 * stepSizeSqr - (3 * node1 - 4 * node2 + node3) * stepSize * x + 2 * (node1 - 2 * node2 + node3) * x * x);
#if 0
					cout << "Iteration " << it << ":  a=" << a << ", b=" << b
						<< ", x=" << x << ", fx=" << fx << ", dfx=" << dfx << endl;
#endif

				x -= fx / dfx;

				if (EXPECT_NOT_TAKEN(x <= a || x >= b || dfx == 0))
					x = 0.5f * (b + a);

				/* Integrated version of the above Lagrange polynomial */
				Float intval = integratedDensity + temp * (1.0f / 6.0f) * (x * (6 * node1 * stepSizeSqr - 3 * (3 * node1 - 4 * node2 + node3) * stepSize * x + 4 * (node1 - 2 * node2 + node3) * x * x));
				fx = intval - desiredDensity;

				if (std::abs(fx) < 1e-6f) {
					t = mint + stepSize * i + x;
					integratedDensity = intval;
					densityAtT = temp * (node1 * stepSizeSqr - (3 * node1 - 4 * node2 + node3) * stepSize * x + 2 * (node1 - 2 * node2 + node3) * x * x);
					return true;
				} else if (++it > 30) {
					Log(EWarn, "invertDensityIntegral(): stuck in Newton-Bisection -- "
							   "round-off error issues? The step size was %e, fx=%f, dfx=%f, "
							   "a=%f, b=%f",
						stepSize, fx, dfx, a, b);
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
			Log(EWarn, "invertDensityIntegral(): unable to make forward progress -- "
					   "round-off error issues? The step size was %e",
				stepSize);
			break;
		}
		integratedDensity = newDensity;
		node1 = node3;
		p = next;
	}

	return false;
}

Spectrum HairMedium::evalTransmittance(const Ray &ray, Sampler *sampler) const {
	if (m_method == ESimpsonQuadrature || sampler == NULL) {
		return Spectrum(math::fastexp(-integrateDensity(ray)));
	} else {
		/* When Woodcock tracking is selected as the sampling method,
		   we can use this method to get a noisy (but unbiased) estimate
		   of the transmittance */
		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return Spectrum(1.0f);
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);

#if defined(HETVOL_STATISTICS)
		avgRayMarchingStepsTransmittance.incrementBase();
#endif
		int nSamples = 2; /// XXX make configurable
		Float result = 0;

		for (int i = 0; i < nSamples; ++i) {
			Float t = mint;
			while (true) {
				t -= math::fastlog(1 - sampler->next1D()) * m_invMaxDensity;
				if (t >= maxt) {
					result += 1;
					break;
				}

				Point p = ray(t);
				Float density = lookupDensity(p, ray.d) * m_scale;

#if defined(HETVOL_STATISTICS)
				++avgRayMarchingStepsTransmittance;
#endif

				if (density * m_invMaxDensity > sampler->next1D())
					break;
			}
		}
		return Spectrum(result / nSamples);
	}
}

bool HairMedium::sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
								Sampler *sampler) const {
	Float integratedDensity, densityAtMinT, densityAtT;
	bool success = false;

	if (m_method == ESimpsonQuadrature) {
		Float desiredDensity = -math::fastlog(1 - sampler->next1D());
		if (invertDensityIntegral(ray, desiredDensity, integratedDensity,
								  mRec.t, densityAtMinT, densityAtT)) {
			mRec.p = ray(mRec.t);
			success = true;
			// Spectrum albedo = m_albedo->lookupSpectrum(mRec.p);
			Spectrum albedo = lookupAlbedo(mRec.p, ray.d);
			mRec.sigmaS = albedo * densityAtT;
			mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
			mRec.orientation = m_orientation != NULL
								   ? m_orientation->lookupVector(mRec.p)
								   : Vector(0.0f);
			mRec.variance = m_variance != nullptr ? m_variance->lookupFloat(mRec.p) : Float(0.f);
		}

		Float expVal = math::fastexp(-integratedDensity);
		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * densityAtT;
		mRec.pdfSuccessRev = expVal * densityAtMinT;
		mRec.transmittance = Spectrum(expVal);
		mRec.time = ray.time;
	} else {
		/* The following information is invalid when
		   using Woodcock-tracking */
		mRec.pdfFailure = 1.0f;
		mRec.pdfSuccess = 1.0f;
		mRec.pdfSuccessRev = 1.0f;
		mRec.transmittance = Spectrum(1.0f);
		mRec.time = ray.time;

#if defined(HETVOL_STATISTICS)
		avgRayMarchingStepsSampling.incrementBase();
#endif

		Float mint, maxt;
		if (!m_densityAABB.rayIntersect(ray, mint, maxt))
			return false;
		mint = std::max(mint, ray.mint);
		maxt = std::min(maxt, ray.maxt);

		Float t = mint, densityAtT = 0;
		while (true) {
			t -= math::fastlog(1 - sampler->next1D()) * m_invMaxDensity;
			if (t >= maxt)
				break;

			Point p = ray(t);
			densityAtT = lookupDensity(p, ray.d) * m_scale;
#if defined(HETVOL_STATISTICS)
			++avgRayMarchingStepsSampling;
#endif
			if (densityAtT > 0) {
				if (densityAtT * m_invMaxDensity > sampler->next1D()) {
					mRec.t = t;
					mRec.p = p;
					// Spectrum albedo = m_albedo->lookupSpectrum(p);
					Spectrum albedo = lookupAlbedo(p, ray.d);
					mRec.sigmaS = albedo * densityAtT;
					mRec.sigmaA = Spectrum(densityAtT) - mRec.sigmaS;
					mRec.transmittance = Spectrum(densityAtT != 0.0f ? 1.0f / densityAtT : 0);
					if (!std::isfinite(mRec.transmittance[0])) // prevent rare overflow warnings
						mRec.transmittance = Spectrum(0.0f);
					mRec.orientation = m_orientation != NULL
										   ? m_orientation->lookupVector(p)
										   : Vector(0.0f);
					mRec.variance = m_variance != nullptr ? m_variance->lookupFloat(mRec.p) : Float(0.f);
					mRec.medium = this;
					success = true;
					break;
				}
			}
		}
	}
	mRec.medium = this;

	return success && mRec.pdfSuccess > 0;
}

void HairMedium::eval(const Ray &ray, MediumSamplingRecord &mRec) const {
	if (m_method == ESimpsonQuadrature) {
		Float expVal = math::fastexp(-integrateDensity(ray));
		Float mintDensity = lookupDensity(ray(ray.mint), ray.d) * m_scale;
		Float maxtDensity = 0.0f;
		Spectrum maxtAlbedo(0.0f);
		if (ray.maxt < std::numeric_limits<Float>::infinity()) {
			Point p = ray(ray.maxt);
			maxtDensity = lookupDensity(p, ray.d) * m_scale;
			// maxtAlbedo = m_albedo->lookupSpectrum(p);
			maxtAlbedo = lookupAlbedo(p, ray.d);
		}
		mRec.transmittance = Spectrum(expVal);
		mRec.pdfFailure = expVal;
		mRec.pdfSuccess = expVal * maxtDensity;
		mRec.pdfSuccessRev = expVal * mintDensity;
		mRec.sigmaS = maxtAlbedo * maxtDensity;
		mRec.sigmaA = Spectrum(maxtDensity) - mRec.sigmaS;
		mRec.time = ray.time;
		mRec.medium = this;
	} else {
		Log(EError, "eval(): unsupported integration method!");
	}
}

bool HairMedium::isHomogeneous() const {
	return false;
}

std::string HairMedium::toString() const {
	std::ostringstream oss;
	oss << "HairMedium[" << endl
		<< "  density = " << indent(m_density.toString()) << "," << endl
		<< "  albedo = " << indent(m_albedo.toString()) << "," << endl
		<< "  orientation = " << indent(m_orientation.toString()) << "," << endl
		<< "  stepSize = " << m_stepSize << "," << endl
		<< "  scale = " << m_scale << endl
		<< "]";
	return oss.str();
}

void HairMedium::evalCoefficientsWithDirection(const Point &p, const Vector3 &d, Spectrum &sigmaS, Spectrum &sigmaT, Vector3 &orientation) const {
	Float density = m_density->lookupFloat(p) * m_scale;
	if (m_anisotropicMedium && density != 0) {
		orientation = m_orientation->lookupVector(p);
	}

	if (m_anisotropicMedium && density != 0) {
		if (!orientation.isZero()) {
			density *= m_phaseFunction->sigmaDir(dot(d, orientation));
		}
	}

	if (density < Epsilon) {
		sigmaS = Spectrum(0.0f);
		sigmaT = Spectrum(0.0f);
		return;
	}

	Spectrum albedo = m_albedo->lookupSpectrum(p);
	if (m_anisotropicAlbedo && !albedo.isZero()) {
		if (!orientation.isZero()) {
			albedo *= m_phaseFunction->albedoDir(d, orientation);
		}
	}

	sigmaS = albedo * density;
	sigmaT = Spectrum(density);
}

void HairMedium::evalCoefficients(const Point &p, Spectrum &sigmaS, Spectrum &sigmaT) const {
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

void HairMedium::evalDensityWithDirection(const Point &p, const Vector3 &d, Float &sigmaT) const {
	Float density = m_density->lookupFloat(p) * m_scale;
	Vector3 orientation;
	if (m_anisotropicMedium && density != 0) {
		orientation = m_orientation->lookupVector(p);
	}

	if (m_anisotropicMedium && density != 0) {
		if (!orientation.isZero()) {
			density *= m_phaseFunction->sigmaDir(dot(d, orientation));
		}
	}

	sigmaT = density;
}

Vector3 HairMedium::getOrientation(const Point &p) const {
	if (m_orientation == nullptr) {
		Log(EError, "m_orientation is null!");
		return Vector3(0.f);
	}
	return m_orientation->lookupVector(p);
}

Float HairMedium::getOrientationVariance(const Point &p) const {
	if (m_variance == nullptr) {
		Log(EError, "m_variance is null!");
		return Float(0.f);
	}
	return m_variance->lookupFloat(p);
}

MTS_IMPLEMENT_CLASS_S(HairMedium, false, Medium)
MTS_EXPORT_PLUGIN(HairMedium, "Heterogeneous medium");
MTS_NAMESPACE_END