#pragma once
#if !defined(__HAIRVOXELIZER2_H_)
#define __HAIRVOXELIZER2_H_

#include <cstdint>
#include <fstream>
#include <mitsuba/core/voxelgrid.h>
#include <mitsuba/render/precomputer.h>
#include <src/shapes/bsplinehair.h>

MTS_NAMESPACE_BEGIN

#ifdef _WIN32
#ifdef HairVoxelizer2_EXPORTS
#define HairVoxelizer2API __declspec(dllexport)
#else
#define HairVoxelizer2API __declspec(dllimport)
#endif
#else
#define HairVoxelizer2API __attribute__((visibility("default")))
#endif

struct HairVoxelizer2API HairVoxelizer2Cell {
	Vector3 averageDirection;
	Float variance;
	Float perpAttenuation;
	Vector3 tensorDiag;
	Vector3 tensorOffDiag;

	HairVoxelizer2Cell()
		: averageDirection(0, 0, 0), variance(0), perpAttenuation(0),
		  tensorDiag(0, 0, 0), tensorOffDiag(0, 0, 0) {
	}
};

struct HairVoxelizer2QuantizedDirection {
	uint8_t theta;
	uint8_t phi;
};

class HairVoxelizer2API HairVoxelizer2 : public Precomputer {
  public:
	using HairVoxelGrid = VoxelGrid<HairVoxelizer2Cell, EStorageOrder::RowMajor>;

	/// Construct from scene properties (voxel resolution, sampling step).
	HairVoxelizer2(const Properties &props);
	/// Stream constructor used by Mitsuba's serialization framework.
	HairVoxelizer2(Stream *stream, InstanceManager *manager);

	/// Entry point: finds hair shape and launches DDA voxelization.
	void runPreBuild(Scene *scene) override;

	/// Access final voxel statistics (direction, variance, density).
	inline const HairVoxelGrid &getVoxelData() const {
		return voxelData;
	}

	/// Access world-space voxel side lengths.
	inline const Vector &getVoxelSize() const {
		return voxelSize;
	}

	/// Access the cubic world-space grid bounds.
	inline AABB getGridAABB() const {
		return gridAABB;
	}

	/// Access voxel resolution.
	inline Vector3i getVoxelCount() const {
		return Vector3i(voxelCountX, voxelCountY, voxelCountZ);
	}

	/// Convert a world-space point to clamped voxel index.
	Vector3i getVoxelIndex(const Point &p) const;
	/// Convert a voxel index to its world-space center.
	Point getVoxelCenter(int x, int y, int z) const;

	MTS_DECLARE_CLASS()
  protected:
	/// Core builder: segment sampling + DDA traversal + tensor finalization.
	void buildVoxelDataDDA(ref<BSplineHairShape> bSplineHairShape);
	/// Convert direction vector into Mitsuba qdir 2-byte representation.
	HairVoxelizer2QuantizedDirection quantizeDirection(const Vector3 &dir) const;
	/// Write Mitsuba VOL header.
	void writeVolHeader(std::ofstream &stream, int32_t encodingType, int32_t channels,
						const Vector3i &resolution, const AABB &aabb) const;
	/// Serialize direction / variance / attenuation volume files.
	void serializeHairVoxelGrid(const std::string &hairName, const fs::path &sceneFilePath) const;

	int voxelCountX, voxelCountY, voxelCountZ;
	Float sampleDeltaAlongHair;
	bool useInputSegments;
	AABB gridAABB;
	Vector voxelSize;
	HairVoxelGrid voxelData;
};

MTS_NAMESPACE_END

#endif
