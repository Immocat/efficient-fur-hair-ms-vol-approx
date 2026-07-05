#pragma once
#if !defined(__HAIRVOXELIZER_H_)
#define __HAIRVOXELIZER_H_

#include <cstdint>
#include <fstream>
#include <mitsuba/core/voxelgrid.h>
#include <mitsuba/render/precomputer.h>
#include <src/shapes/bsplinehair.h>

MTS_NAMESPACE_BEGIN

#ifdef _WIN32
#ifdef HairVoxel_EXPORTS
#define HairVoxelAPI __declspec(dllexport)
#else
#define HairVoxelAPI __declspec(dllimport)
#endif
#else
#define HairVoxelAPI __attribute__((visibility("default")))
#endif

struct HairVoxelAPI HairVoxelT {
	Vector3 averageDirection;
	Float variance;
	Float perpAttenuation;

	HairVoxelT()
		: averageDirection(0, 0, 0), variance(0), perpAttenuation(0) {
	}
};

struct HairSamplePoint {
	Point p;
	Vector tangent;
	uint32_t strandID;
};

struct QuantizedDirection {
	uint8_t theta;
	uint8_t phi;
};

/**
 * Precomputes voxelized hair statistics for volume-based hair rendering.
 *
 * The voxelizer samples the input hair strands, gathers nearby samples per
 * occupied voxel, and stores tangent direction, tangent variance, and transverse
 * attenuation for the hair volume model.
 */
class HairVoxelAPI HairVoxel : public Precomputer {
  public:
	using HairVoxelGrid = VoxelGrid<HairVoxelT, EStorageOrder::RowMajor>;

	HairVoxel(const Properties &props);
	HairVoxel(Stream *stream, InstanceManager *manager);

	void runPreBuild(Scene *scene) override;

	inline const HairVoxelGrid &getVoxelData() const {
		return voxelData;
	}

	inline const Vector &getVoxelSize() const {
		return voxelSize;
	}

	inline AABB getGridAABB() const {
		return gridAABB;
	}

	inline Vector3i getVoxelCount() const {
		return Vector3i(voxelCountX, voxelCountY, voxelCountZ);
	}

	Vector3i getVoxelIndex(const Point &p) const;
	Point getVoxelCenter(int x, int y, int z) const;

	MTS_DECLARE_CLASS()

  protected:
	void buildVoxelData(ref<BSplineHairShape> bSplineHairShape);
	QuantizedDirection quantizeDirection(const Vector3 &dir) const;
	void writeVolHeader(std::ofstream &stream, int32_t encodingType, int32_t channels,
						const Vector3i &resolution, const AABB &aabb) const;
	void serializeHairVoxelGrid(const std::string &hairName, const fs::path &sceneFilePath) const;

	int voxelCountX, voxelCountY, voxelCountZ;
	Float dMultiplier;
	Float sampleDeltaAlongHair;
	std::string outputSuffix;
	AABB gridAABB;
	Vector voxelSize;
	HairVoxelGrid voxelData;
};

MTS_NAMESPACE_END
#endif
