#pragma once
#include <cstdint>
#if !defined(__MITSUBA_RENDER_RAYMASK_H_)
#define __MITSUBA_RENDER_RAYMASK_H_

using RayMaskType = uint32_t;

enum ERayMaskCommon : RayMaskType {
	ENone = 0,
	EDefault = 1,
};

enum EHairRayMask : RayMaskType {
	EHairGeometry = 2,
	EHairVolumeRed = 4,
	EHairVolumeGreen = 8,
	EHairVolumeBlue = 16,
	EHairVolume = EHairVolumeRed | EHairVolumeGreen | EHairVolumeBlue,

};

inline bool masked(RayMaskType rayMask, RayMaskType shapeRayMask) {
	return rayMask & shapeRayMask;
}

inline bool notMasked(RayMaskType rayMask, RayMaskType shapeRayMask) {
	return !masked(rayMask, shapeRayMask);
}

#endif