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

#ifndef __MITSUBA_RENDER_HAIRMEDIUM_H_
#define __MITSUBA_RENDER_HAIRMEDIUM_H_

#include <mitsuba/core/shvector.h>
#include <mitsuba/core/voxelgrid.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/volume.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Heterogeneous participating medium implementation
 * (Referred to as HairMedium in the codebase)
 */
class HairMedium : public Medium {
  public:
	/// Possible integration modes
	enum EIntegrationMethod {
		/**
		 * \brief Use deterministic composite Simpson quadrature both
		 * to compute transmittances, and to sample scattering locations
		 */
		ESimpsonQuadrature = 0,

		/**
		 * \brief Use stochastic Woodcock tracking. This is potentially
		 * faster and more robust, but has the disadvantage of being
		 * incompatible with bidirectional rendering methods, which
		 * usually need to know the probability of a sample.
		 */
		EWoodcockTracking
	};

	HairMedium(const Properties &props);
	HairMedium(Stream *stream, InstanceManager *manager);

	virtual ~HairMedium();

	/* Serialize the volume to a binary data stream */
	void serialize(Stream *stream, InstanceManager *manager) const;

	void configure();

	void addChild(const std::string &name, ConfigurableObject *child);

	/*
	 * This function uses Simpson quadrature to compute following
	 * integral:
	 *
	 *    \int_{ray.mint}^{ray.maxt} density(ray(x)) dx
	 */
	Float integrateDensity(const Ray &ray) const;

	/**
	 * This function uses composite Simpson quadrature to solve the
	 * following integral equation for \a t:
	 *
	 *    \int_{ray.mint}^t density(ray(x)) dx == desiredDensity
	 */
	bool invertDensityIntegral(const Ray &ray, Float desiredDensity,
							   Float &integratedDensity, Float &t, Float &densityAtMinT,
							   Float &densityAtT) const;

	Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const;

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
						Sampler *sampler) const;

	void eval(const Ray &ray, MediumSamplingRecord &mRec) const;

	bool isHomogeneous() const;

	std::string toString() const;

	void evalCoefficientsWithDirection(const Point &p, const Vector3 &d, Spectrum &sigmaS, Spectrum &sigmaT, Vector3 &orientation) const override;

	void evalCoefficients(const Point &p, Spectrum &sigmaS, Spectrum &sigmaT) const override;

	void evalDensityWithDirection(const Point &p, const Vector3 &d, Float &sigmaT) const override;

	MTS_DECLARE_CLASS()

  protected:
	inline Float lookupDensity(const Point &p, const Vector &d) const;

	inline Spectrum lookupAlbedo(const Point &p, const Vector &d) const;

	Vector3 getOrientation(const Point &p) const override;

	Float getOrientationVariance(const Point &p) const override;

  protected:
	EIntegrationMethod m_method;
	ref<VolumeDataSource> m_density;
	ref<VolumeDataSource> m_albedo;
	ref<VolumeDataSource> m_orientation;
	ref<VolumeDataSource> m_variance;

	VoxelGrid<SHVectorRGB> m_LPV;
	int m_LPVRes;

	Float m_scale;
	bool m_anisotropicMedium;
	bool m_anisotropicAlbedo;
	Float m_stepSize;
	AABB m_densityAABB;
	Float m_maxDensity;
	Float m_invMaxDensity;

	bool m_biasedStepping;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_HAIRMEDIUM_H_ */