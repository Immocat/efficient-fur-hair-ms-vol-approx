#include "hairvoxelizer2.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

namespace {

/// Per-voxel temporary accumulators used before tensor finalization.
struct VoxelAccum {
	Float weight;
	Matrix3x3 moment;

	VoxelAccum()
		: weight(0.0f), moment(0.0f) {
	}
};

/// Convert 3D voxel index to linear storage index.
inline uint32_t flattenVoxelIndex(int x, int y, int z, int nx, int ny) {
	return (uint32_t)(x + nx * (y + ny * z));
}

/// Convert linear storage index back to 3D voxel index.
inline Vector3i unflattenVoxelIndex(uint32_t idx, int nx, int ny) {
	int z = (int)(idx / (uint32_t)(nx * ny));
	uint32_t rem = idx - (uint32_t)(z * nx * ny);
	int y = (int)(rem / (uint32_t)nx);
	int x = (int)(rem - (uint32_t)(y * nx));
	return Vector3i(x, y, z);
}

/// Clamp scalar to [0, 1].
inline Float clamp01(Float x) {
	return std::max((Float)0.0f, std::min((Float)1.0f, x));
}

inline Float lerpScalar(Float a, Float b, Float t) {
	return a + (b - a) * t;
}

inline Vector mulSymmetricTensor(Float txx, Float tyy, Float tzz,
								 Float txy, Float txz, Float tyz, const Vector &v) {
	return Vector(
		txx * v.x + txy * v.y + txz * v.z,
		txy * v.x + tyy * v.y + tyz * v.z,
		txz * v.x + tyz * v.y + tzz * v.z);
}

/// Falcor parity: dominant direction by power iteration from initial (1,0,0).
inline Vector principalDirectionFromTensor(Float txx, Float tyy, Float tzz,
										   Float txy, Float txz, Float tyz) {
	Vector v(1.0f, 0.0f, 0.0f);
	for (int i = 0; i < 8; ++i) {
		Vector Tv = mulSymmetricTensor(txx, tyy, tzz, txy, txz, tyz, v);
		Float len = Tv.length();
		v = (len > 1e-6f) ? (Tv / len) : Vector(1.0f, 0.0f, 0.0f);
	}
	Float len2 = v.lengthSquared();
	return (len2 > 1e-12f) ? normalize(v) : Vector(0.0f);
}

/// Liang-Barsky style segment/AABB clipping in parametric t \in [0, 1].
bool clipSegmentToAABB(const Point &a, const Point &b, const AABB &box, Float &tMin, Float &tMax) {
	tMin = 0.0f;
	tMax = 1.0f;
	Vector d = b - a;

	for (int axis = 0; axis < 3; ++axis) {
		Float p0 = a[axis];
		Float dir = d[axis];
		Float minv = box.min[axis];
		Float maxv = box.max[axis];

		if (std::abs(dir) < 1e-12f) {
			if (p0 < minv || p0 > maxv)
				return false;
			continue;
		}

		Float invDir = 1.0f / dir;
		Float t1 = (minv - p0) * invDir;
		Float t2 = (maxv - p0) * invDir;
		if (t1 > t2)
			std::swap(t1, t2);

		tMin = std::max(tMin, t1);
		tMax = std::min(tMax, t2);

		if (tMin > tMax)
			return false;
	}

	return true;
}

/// Traverse one clipped segment through crossed voxels and accumulate
/// - projected-area weight W contribution
/// - structure tensor moments M = sum(w * v v^T).
void accumulateSegmentDDA(const Point &segA,
						  const Point &segB,
						  Float radiusA,
						  Float radiusB,
						  const AABB &gridAABB,
						  const Vector &voxelSize,
						  int voxelCountX,
						  int voxelCountY,
						  int voxelCountZ,
						  std::unordered_map<uint32_t, VoxelAccum> &accum) {
	Vector segment = segB - segA;
	Float segmentLen = segment.length();
	if (segmentLen <= 1e-10f)
		return;

	Float clipTMin, clipTMax;
	if (!clipSegmentToAABB(segA, segB, gridAABB, clipTMin, clipTMax))
		return;

	clipTMin = clamp01(clipTMin);
	clipTMax = clamp01(clipTMax);
	if (clipTMax <= clipTMin)
		return;

	Point pIn = segA + segment * clipTMin;
	Point pOut = segA + segment * clipTMax;
	Float rIn = lerpScalar(radiusA, radiusB, clipTMin);
	Float rOut = lerpScalar(radiusA, radiusB, clipTMax);

	Vector dirWorld = pOut - pIn;
	Float segmentLengthWorld = dirWorld.length();
	if (segmentLengthWorld <= 1e-10f)
		return;
	Vector direction = dirWorld / segmentLengthWorld;

	Vector pInGrid = CwiseDiv((pIn - gridAABB.min), voxelSize);
	Vector pOutGrid = CwiseDiv((pOut - gridAABB.min), voxelSize);
	Vector dirGrid = pOutGrid - pInGrid;

	int ix = math::clamp(math::floorToInt(pInGrid.x), 0, voxelCountX - 1);
	int iy = math::clamp(math::floorToInt(pInGrid.y), 0, voxelCountY - 1);
	int iz = math::clamp(math::floorToInt(pInGrid.z), 0, voxelCountZ - 1);

	int stepX = 0, stepY = 0, stepZ = 0;
	Float tMaxX = std::numeric_limits<Float>::infinity();
	Float tMaxY = std::numeric_limits<Float>::infinity();
	Float tMaxZ = std::numeric_limits<Float>::infinity();
	Float tDeltaX = std::numeric_limits<Float>::infinity();
	Float tDeltaY = std::numeric_limits<Float>::infinity();
	Float tDeltaZ = std::numeric_limits<Float>::infinity();

	if (dirGrid.x > 1e-12f) {
		stepX = 1;
		Float nextBoundary = (Float)ix + 1.0f;
		tMaxX = (nextBoundary - pInGrid.x) / dirGrid.x;
		tDeltaX = std::abs(1.0f / dirGrid.x);
	} else if (dirGrid.x < -1e-12f) {
		stepX = -1;
		Float nextBoundary = (Float)ix;
		tMaxX = (nextBoundary - pInGrid.x) / dirGrid.x;
		tDeltaX = std::abs(1.0f / dirGrid.x);
	}

	if (dirGrid.y > 1e-12f) {
		stepY = 1;
		Float nextBoundary = (Float)iy + 1.0f;
		tMaxY = (nextBoundary - pInGrid.y) / dirGrid.y;
		tDeltaY = std::abs(1.0f / dirGrid.y);
	} else if (dirGrid.y < -1e-12f) {
		stepY = -1;
		Float nextBoundary = (Float)iy;
		tMaxY = (nextBoundary - pInGrid.y) / dirGrid.y;
		tDeltaY = std::abs(1.0f / dirGrid.y);
	}

	if (dirGrid.z > 1e-12f) {
		stepZ = 1;
		Float nextBoundary = (Float)iz + 1.0f;
		tMaxZ = (nextBoundary - pInGrid.z) / dirGrid.z;
		tDeltaZ = std::abs(1.0f / dirGrid.z);
	} else if (dirGrid.z < -1e-12f) {
		stepZ = -1;
		Float nextBoundary = (Float)iz;
		tMaxZ = (nextBoundary - pInGrid.z) / dirGrid.z;
		tDeltaZ = std::abs(1.0f / dirGrid.z);
	}

	Float tCurrent = 0.0f;
	const Float tEnd = 1.0f;
	const Float tEps = 1e-7f;
	const int maxSteps = voxelCountX + voxelCountY + voxelCountZ + 8;
	int iter = 0;

	while (ix >= 0 && ix < voxelCountX &&
		   iy >= 0 && iy < voxelCountY &&
		   iz >= 0 && iz < voxelCountZ &&
		   tCurrent < tEnd && iter++ < maxSteps) {

		Float tMin = std::min(tMaxX, std::min(tMaxY, tMaxZ));
		Float tNext = std::min(tEnd, tMin);
		Float dt = tNext - tCurrent;

		if (dt > tEps) {
			uint32_t linearIdx = flattenVoxelIndex(ix, iy, iz, voxelCountX, voxelCountY);
			VoxelAccum &entry = accum[linearIdx];

			Float tMid = 0.5f * (tCurrent + tNext);
			Float radius = std::max((Float)0.0f, lerpScalar(rIn, rOut, tMid));
			Float weight = dt * segmentLengthWorld * (2.0f * radius);
			if (weight > 0.0f) {
				entry.weight += weight;
				entry.moment(0, 0) += weight * direction.x * direction.x;
				entry.moment(0, 1) += weight * direction.x * direction.y;
				entry.moment(0, 2) += weight * direction.x * direction.z;
				entry.moment(1, 0) += weight * direction.y * direction.x;
				entry.moment(1, 1) += weight * direction.y * direction.y;
				entry.moment(1, 2) += weight * direction.y * direction.z;
				entry.moment(2, 0) += weight * direction.z * direction.x;
				entry.moment(2, 1) += weight * direction.z * direction.y;
				entry.moment(2, 2) += weight * direction.z * direction.z;
			}
		}

		if (tMin > 1.0f)
			break;

		if (tMaxX <= tMin + tEps) {
			ix += stepX;
			tMaxX += tDeltaX;
		}
		if (tMaxY <= tMin + tEps) {
			iy += stepY;
			tMaxY += tDeltaY;
		}
		if (tMaxZ <= tMin + tEps) {
			iz += stepZ;
			tMaxZ += tDeltaZ;
		}

		tCurrent = tNext;
	}
}

} // namespace

/// Parse voxelization parameters from plugin properties.
HairVoxelizer2::HairVoxelizer2(const Properties &props)
	: Precomputer(props) {
	int voxelCount = props.getInteger("voxelCount", 128);
	voxelCountX = props.getInteger("voxelCountX", voxelCount);
	voxelCountY = props.getInteger("voxelCountY", voxelCount);
	voxelCountZ = props.getInteger("voxelCountZ", voxelCount);

	if (voxelCountX == 0)
		voxelCountX = voxelCount;
	if (voxelCountY == 0)
		voxelCountY = voxelCount;
	if (voxelCountZ == 0)
		voxelCountZ = voxelCount;

	sampleDeltaAlongHair = props.getFloat("sampleDeltaAlongHair", 0.02f);
	if (sampleDeltaAlongHair <= 0.0f) {
		Log(EWarn, "HairVoxelizer2: sampleDeltaAlongHair must be > 0. Resetting to 0.02.");
		sampleDeltaAlongHair = 0.02f;
	}
	useInputSegments = props.getBoolean("useInputSegments", false);

	voxelData = HairVoxelGrid(voxelCountX, voxelCountY, voxelCountZ);
}

/// Serialization constructor (kept minimal for this v2 plugin).
HairVoxelizer2::HairVoxelizer2(Stream *stream, InstanceManager *manager)
	: Precomputer(stream, manager) {
}

/// Locate the target hair shape and run the DDA voxel precompute.
void HairVoxelizer2::runPreBuild(Scene *scene) {
	Log(EInfo, "Precomputing HairVoxelizer2.");
	ref<Timer> timer = new Timer();

	ref<BSplineHairShape> bsplinehairshapeRef;
	auto shapes = scene->getShapes();
	for (auto &shape : shapes) {
		BSplineHairShape *bsplinehairshape = dynamic_cast<BSplineHairShape *>(shape.get());
		if (!bsplinehairshape)
			continue;

		if (bsplinehairshapeRef) {
			Log(EError, "HairVoxelizer2 supports exactly one bSpline hair shape.");
			return;
		}

		bsplinehairshapeRef = ref<BSplineHairShape>(bsplinehairshape);
	}

	if (!bsplinehairshapeRef) {
		Log(EError, "HairVoxelizer2 requires a bSpline hair shape in the scene.");
		return;
	}

	buildVoxelDataDDA(bsplinehairshapeRef);

	std::string hairFileName = bsplinehairshapeRef->getFileName();
	size_t dotPos = hairFileName.find_last_of('.');
	std::string hairName = (dotPos != std::string::npos) ? hairFileName.substr(0, dotPos) : hairFileName;
	serializeHairVoxelGrid(hairName, scene->getSourceFile());

	Log(EInfo, "Precomputing HairVoxelizer2 done. Time taken: %i ms.", timer->getMilliseconds());
}

/// World point -> clamped voxel index conversion.
Vector3i HairVoxelizer2::getVoxelIndex(const Point &p) const {
	Vector3 rel = CwiseDiv((p - gridAABB.min), voxelSize);
	return Vector3i(
		math::clamp(math::floorToInt(rel.x), 0, voxelCountX - 1),
		math::clamp(math::floorToInt(rel.y), 0, voxelCountY - 1),
		math::clamp(math::floorToInt(rel.z), 0, voxelCountZ - 1));
}

/// Voxel index -> world-space voxel center.
Point HairVoxelizer2::getVoxelCenter(int x, int y, int z) const {
	return gridAABB.min + CwiseMul(voxelSize, Vector((Float)x + 0.5f, (Float)y + 0.5f, (Float)z + 0.5f));
}

/// Main algorithm:
/// 1) build cubic grid around hair AABB
/// 2) sample curve into line segments
/// 3) DDA-accumulate per-thread sparse voxel stats
/// 4) merge thread-local accumulators
/// 5) finalize average direction, variance, and density per voxel.
void HairVoxelizer2::buildVoxelDataDDA(ref<BSplineHairShape> bSplineHairShape) {
	size_t numPrimitives = bSplineHairShape->getPrimitiveCount();
	if (numPrimitives == 0 && !useInputSegments) {
		Log(EError, "HairVoxelizer2: the hair shape contains no primitives.");
		return;
	}

	const std::vector<Point> &inputVertices = bSplineHairShape->getVertices();
	const std::vector<bool> &inputStartsFiber = bSplineHairShape->getStartFiber();
	const Float inputRadius = bSplineHairShape->getRadius();
	std::vector<size_t> inputSegmentStarts;
	size_t workCount = numPrimitives;

	if (useInputSegments) {
		if (inputVertices.size() < 2 || inputStartsFiber.size() < inputVertices.size()) {
			Log(EError, "HairVoxelizer2: invalid input strand layout for useInputSegments=true.");
			return;
		}
		inputSegmentStarts.reserve(inputVertices.size());
		for (size_t i = 0; i + 1 < inputVertices.size(); ++i) {
			if (!inputStartsFiber[i + 1])
				inputSegmentStarts.push_back(i);
		}
		workCount = inputSegmentStarts.size();
		if (workCount == 0) {
			Log(EError, "HairVoxelizer2: no valid input segments found for useInputSegments=true.");
			return;
		}
	}

	AABB tightAABB = bSplineHairShape->getAABB();
	Vector extents = tightAABB.getExtents();
	Point center = tightAABB.getCenter();
	Float maxSide = std::max(extents.x, std::max(extents.y, extents.z));
	Vector halfSideVector(maxSide * 0.5f);
	gridAABB = AABB(center - halfSideVector, center + halfSideVector);

	voxelSize = Vector(
		gridAABB.getExtents().x / (Float)voxelCountX,
		gridAABB.getExtents().y / (Float)voxelCountY,
		gridAABB.getExtents().z / (Float)voxelCountZ);

	Float voxelVolume = voxelSize.x * voxelSize.y * voxelSize.z;
	if (voxelVolume <= 0.0f) {
		Log(EError, "HairVoxelizer2: invalid voxel volume.");
		return;
	}

	unsigned int hwThreads = std::max(1u, std::thread::hardware_concurrency());
	size_t threadCount = std::min((size_t)hwThreads, workCount);
	threadCount = std::max((size_t)1, threadCount);

	std::vector<std::unordered_map<uint32_t, VoxelAccum>> localVoxelAccums(threadCount);
	std::vector<std::thread> workers;
	workers.reserve(threadCount);

	if (useInputSegments)
		Log(EInfo, "HairVoxelizer2: DDA traversal over %zu input segments using %zu threads.", workCount, threadCount);
	else
		Log(EInfo, "HairVoxelizer2: DDA traversal over %zu primitives using %zu threads.", workCount, threadCount);

	const Float du = sampleDeltaAlongHair;
	for (size_t tid = 0; tid < threadCount; ++tid) {
		size_t start = (workCount * tid) / threadCount;
		size_t end = (workCount * (tid + 1)) / threadCount;

		workers.emplace_back([&, tid, start, end]() {
			auto &localAccum = localVoxelAccums[tid];
			localAccum.reserve(std::max((size_t)1024, (end - start) * 8));

			// Each thread processes a disjoint primitive range into a local sparse map.
			if (useInputSegments) {
				for (size_t segIdx = start; segIdx < end; ++segIdx) {
					size_t vtx = inputSegmentStarts[segIdx];
					accumulateSegmentDDA(
						inputVertices[vtx], inputVertices[vtx + 1],
						inputRadius, inputRadius,
						gridAABB, voxelSize,
						voxelCountX, voxelCountY, voxelCountZ,
						localAccum);
				}
			} else {
				for (size_t primID = start; primID < end; ++primID) {
					CurvePointInfo prevInfo;
					if (!bSplineHairShape->evaluatePrimitive((uint32_t)primID, 0.0f, prevInfo))
						continue;

					Float u = 0.0f;
					while (u < 1.0f) {
						Float nextU = std::min((Float)1.0f, u + du);
						CurvePointInfo nextInfo;
						if (!bSplineHairShape->evaluatePrimitive((uint32_t)primID, nextU, nextInfo))
							break;

						accumulateSegmentDDA(
							prevInfo.p, nextInfo.p,
							prevInfo.radius, nextInfo.radius,
							gridAABB, voxelSize,
							voxelCountX, voxelCountY, voxelCountZ,
							localAccum);

						prevInfo = nextInfo;
						u = nextU;
					}
				}
			}
		});
	}

	for (auto &worker : workers)
		worker.join();

	size_t voxelTotal = (size_t)voxelCountX * voxelCountY * voxelCountZ;
	std::vector<VoxelAccum> merged(voxelTotal);
	std::vector<bool> touched(voxelTotal, false);

	for (size_t tid = 0; tid < threadCount; ++tid) {
		const auto &localAccum = localVoxelAccums[tid];
		for (const auto &kv : localAccum) {
			uint32_t idx = kv.first;
			const VoxelAccum &src = kv.second;
			VoxelAccum &dst = merged[idx];
			dst.weight += src.weight;
			dst.moment += src.moment;
			touched[idx] = true;
		}
	}

	voxelData = HairVoxelGrid(voxelCountX, voxelCountY, voxelCountZ);
	for (size_t linearIdx = 0; linearIdx < voxelTotal; ++linearIdx) {
		if (!touched[linearIdx])
			continue;

		const VoxelAccum &acc = merged[linearIdx];
		if (acc.weight <= 0.0f)
			continue;

		const Float sigmaBase = acc.weight / voxelVolume;
		if (sigmaBase <= 0.0f)
			continue;

		const Float invW = 1.0f / acc.weight;
		Float txx = acc.moment(0, 0) * invW;
		Float tyy = acc.moment(1, 1) * invW;
		Float tzz = acc.moment(2, 2) * invW;
		Float txy = acc.moment(0, 1) * invW;
		Float txz = acc.moment(0, 2) * invW;
		Float tyz = acc.moment(1, 2) * invW;

		const Float trace = std::max((Float)1e-10f, txx + tyy + tzz);
		const Float invTrace = 1.0f / trace;
		txx *= invTrace;
		tyy *= invTrace;
		tzz *= invTrace;
		txy *= invTrace;
		txz *= invTrace;
		tyz *= invTrace;

		Vector principal = principalDirectionFromTensor(txx, tyy, tzz, txy, txz, tyz);
		Float variance = 0.0f;
		if (!principal.isZero()) {
			const Vector w = normalize(principal);
			const Float pParallel =
				txx * w.x * w.x +
				tyy * w.y * w.y +
				tzz * w.z * w.z +
				2.0f * txy * w.x * w.y +
				2.0f * txz * w.x * w.z +
				2.0f * tyz * w.y * w.z;
			variance = clamp01(1.0f - clamp01(pParallel));
		}

		Vector3i idx = unflattenVoxelIndex((uint32_t)linearIdx, voxelCountX, voxelCountY);
		HairVoxelizer2Cell &cell = voxelData(idx.x, idx.y, idx.z);
		cell.averageDirection = principal;
		cell.variance = variance;
		cell.perpAttenuation = sigmaBase;
		cell.tensorDiag = Vector3(txx, tyy, tzz);
		cell.tensorOffDiag = Vector3(txy, txz, tyz);
	}
}

HairVoxelizer2QuantizedDirection HairVoxelizer2::quantizeDirection(const Vector3 &dir) const {
	if (dir.lengthSquared() < 1e-8f)
		return {255, 0};

	Vector3 n_dir = normalize(dir);
	float clamped_z = std::max(-1.0f, std::min(1.0f, (float)n_dir.z));
	float theta = std::acos(clamped_z);
	float phi = std::atan2((float)n_dir.y, (float)n_dir.x);
	if (phi < 0.0f)
		phi += 2.0f * (float)M_PI;

	uint8_t theta_byte = static_cast<uint8_t>(std::round((theta / (float)M_PI) * 254.0f));
	uint8_t phi_byte = static_cast<uint8_t>(std::round((phi / (2.0f * (float)M_PI)) * 254.0f));
	return {theta_byte, phi_byte};
}

void HairVoxelizer2::writeVolHeader(std::ofstream &stream, int32_t encodingType, int32_t channels,
									const Vector3i &resolution, const AABB &aabb) const {
	stream.write("VOL", 3);

	uint8_t version = 3;
	stream.write(reinterpret_cast<const char *>(&version), 1);

	stream.write(reinterpret_cast<const char *>(&encodingType), sizeof(int32_t));
	stream.write(reinterpret_cast<const char *>(&resolution.x), sizeof(int32_t));
	stream.write(reinterpret_cast<const char *>(&resolution.y), sizeof(int32_t));
	stream.write(reinterpret_cast<const char *>(&resolution.z), sizeof(int32_t));
	stream.write(reinterpret_cast<const char *>(&channels), sizeof(int32_t));

	float xmin = (float)aabb.min.x, ymin = (float)aabb.min.y, zmin = (float)aabb.min.z;
	float xmax = (float)aabb.max.x, ymax = (float)aabb.max.y, zmax = (float)aabb.max.z;
	stream.write(reinterpret_cast<const char *>(&xmin), sizeof(float));
	stream.write(reinterpret_cast<const char *>(&ymin), sizeof(float));
	stream.write(reinterpret_cast<const char *>(&zmin), sizeof(float));
	stream.write(reinterpret_cast<const char *>(&xmax), sizeof(float));
	stream.write(reinterpret_cast<const char *>(&ymax), sizeof(float));
	stream.write(reinterpret_cast<const char *>(&zmax), sizeof(float));
}

void HairVoxelizer2::serializeHairVoxelGrid(const std::string &hairName, const fs::path &sceneFilePath) const {
	const Vector3i &res = voxelData.getResolution();
	if (voxelData.getSize() == 0) {
		Log(EError, "HairVoxelizer2: Attempting to serialize an empty VoxelGrid.");
		return;
	}

	boost::filesystem::path outputDirectory = sceneFilePath.parent_path();
	std::stringstream ss;
	ss << hairName << "_" << res.x << "_" << res.y << "_" << res.z;
	std::string baseFilename = ss.str();

	// 1) Average direction (qdir)
	boost::filesystem::path dirPath = outputDirectory / (baseFilename + "_dir.vol");
	std::ofstream dirStream(dirPath.string(), std::ios::binary);
	if (!dirStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open file for writing: \"%s\"", dirPath.string().c_str());
		return;
	}
	writeVolHeader(dirStream, 4, 3, res, gridAABB);
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x) {
				HairVoxelizer2QuantizedDirection qdir = quantizeDirection(voxelData(x, y, z).averageDirection);
				dirStream.write(reinterpret_cast<const char *>(&qdir.theta), 1);
				dirStream.write(reinterpret_cast<const char *>(&qdir.phi), 1);
			}
		}
	}
	dirStream.close();

	// 2) Variance
	boost::filesystem::path varPath = outputDirectory / (baseFilename + "_var.vol");
	std::ofstream varStream(varPath.string(), std::ios::binary);
	if (!varStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open file for writing: \"%s\"", varPath.string().c_str());
		return;
	}
	writeVolHeader(varStream, 1, 1, res, gridAABB);
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x) {
				float variance = (float)voxelData(x, y, z).variance;
				varStream.write(reinterpret_cast<const char *>(&variance), sizeof(float));
			}
		}
	}
	varStream.close();

	// 3) Tensor diagonal channels: (Txx, Tyy, Tzz)
	boost::filesystem::path tensorDiagPath = outputDirectory / (baseFilename + "_tensorDiag.vol");
	std::ofstream tensorDiagStream(tensorDiagPath.string(), std::ios::binary);
	if (!tensorDiagStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open file for writing: \"%s\"", tensorDiagPath.string().c_str());
		return;
	}
	writeVolHeader(tensorDiagStream, 1, 3, res, gridAABB);
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x) {
				const Vector3 &diag = voxelData(x, y, z).tensorDiag;
				const float txx = (float)diag.x;
				const float tyy = (float)diag.y;
				const float tzz = (float)diag.z;
				tensorDiagStream.write(reinterpret_cast<const char *>(&txx), sizeof(float));
				tensorDiagStream.write(reinterpret_cast<const char *>(&tyy), sizeof(float));
				tensorDiagStream.write(reinterpret_cast<const char *>(&tzz), sizeof(float));
			}
		}
	}
	tensorDiagStream.close();

	// 4) Tensor off-diagonal channels: (Txy, Txz, Tyz)
	boost::filesystem::path tensorOffDiagPath = outputDirectory / (baseFilename + "_tensorOffDiag.vol");
	std::ofstream tensorOffDiagStream(tensorOffDiagPath.string(), std::ios::binary);
	if (!tensorOffDiagStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open file for writing: \"%s\"", tensorOffDiagPath.string().c_str());
		return;
	}
	writeVolHeader(tensorOffDiagStream, 1, 3, res, gridAABB);
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x) {
				const Vector3 &off = voxelData(x, y, z).tensorOffDiag;
				const float txy = (float)off.x;
				const float txz = (float)off.y;
				const float tyz = (float)off.z;
				tensorOffDiagStream.write(reinterpret_cast<const char *>(&txy), sizeof(float));
				tensorOffDiagStream.write(reinterpret_cast<const char *>(&txz), sizeof(float));
				tensorOffDiagStream.write(reinterpret_cast<const char *>(&tyz), sizeof(float));
			}
		}
	}
	tensorOffDiagStream.close();

	// 5) Attenuation (raw .hairvol + normalized .vol) + scale xml
	boost::filesystem::path rawHairVolPath = outputDirectory / (baseFilename + "_atten_raw.hairvol");
	std::ofstream rawHairVolStream(rawHairVolPath.string(), std::ios::binary);
	if (!rawHairVolStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open file for writing: \"%s\"", rawHairVolPath.string().c_str());
		return;
	}
	writeVolHeader(rawHairVolStream, 1, 1, res, gridAABB);
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x) {
				float rawDensity = (float)voxelData(x, y, z).perpAttenuation;
				rawHairVolStream.write(reinterpret_cast<const char *>(&rawDensity), sizeof(float));
			}
		}
	}
	rawHairVolStream.close();
	Log(EInfo, "HairVoxelizer2: Wrote raw attenuation .hairvol to \"%s\"", rawHairVolPath.string().c_str());

	float maxDensity = 0.0f;
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x)
				maxDensity = std::max(maxDensity, (float)voxelData(x, y, z).perpAttenuation);
		}
	}
	if (maxDensity <= 1e-9f) {
		maxDensity = 1.0f;
		Log(EWarn, "HairVoxelizer2: Max density is close to 0, setting scale to 1.0.");
	}

	ss.str("");
	ss.clear();
	ss << std::fixed << std::setprecision(4) << maxDensity;
	std::string rawFloatStr = ss.str();
	if (rawFloatStr.find('.') != std::string::npos) {
		size_t lastNonZero = rawFloatStr.find_last_not_of('0');
		if (lastNonZero != std::string::npos && rawFloatStr[lastNonZero] == '.')
			rawFloatStr.erase(lastNonZero);
		else if (lastNonZero != std::string::npos)
			rawFloatStr.erase(lastNonZero + 1);
	}
	std::string filenameSuffix = rawFloatStr;
	std::replace(filenameSuffix.begin(), filenameSuffix.end(), '.', '-');

	std::string attenFilename = baseFilename + "_atten_" + filenameSuffix + ".vol";
	boost::filesystem::path attenPath = outputDirectory / attenFilename;
	std::ofstream attenStream(attenPath.string(), std::ios::binary);
	if (!attenStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open file for writing: \"%s\"", attenPath.string().c_str());
		return;
	}
	writeVolHeader(attenStream, 1, 1, res, gridAABB);
	for (int z = 0; z < res.z; ++z) {
		for (int y = 0; y < res.y; ++y) {
			for (int x = 0; x < res.x; ++x) {
				float normalizedDensity = (float)voxelData(x, y, z).perpAttenuation / maxDensity;
				normalizedDensity = std::max(0.0f, std::min(1.0f, normalizedDensity));
				attenStream.write(reinterpret_cast<const char *>(&normalizedDensity), sizeof(float));
			}
		}
	}
	attenStream.close();

	boost::filesystem::path xmlPath = outputDirectory / (hairName + "_densityScale.xml");
	std::ofstream xmlStream(xmlPath.string());
	if (!xmlStream.is_open()) {
		Log(EError, "HairVoxelizer2: Could not open XML file for writing: \"%s\"", xmlPath.string().c_str());
		return;
	}
	xmlStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
	xmlStream << "<scene version=\"0.6.0\">" << std::endl;
	xmlStream << "    <default name=\"DensityScale\" value=\"" << filenameSuffix << "\"/>" << std::endl;
	xmlStream << "    <default name=\"DensityValue\" value=\"" << rawFloatStr << "\"/>" << std::endl;
	xmlStream << "</scene>" << std::endl;
	xmlStream.close();

	Log(EInfo, "HairVoxelizer2: Serialized voxel volumes to \"%s\"", outputDirectory.string().c_str());
}

MTS_IMPLEMENT_CLASS_API(HairVoxelizer2, false, Precomputer, HairVoxelizer2API)
MTS_EXPORT_PLUGIN(HairVoxelizer2, "Hair voxelization precomputer v2 (DDA)");

MTS_NAMESPACE_END
