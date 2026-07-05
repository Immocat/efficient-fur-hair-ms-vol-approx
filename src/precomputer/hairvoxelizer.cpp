#include "hairvoxelizer.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <vector>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

using HairSampleVoxelGrid = VoxelGrid<std::vector<HairSamplePoint>, EStorageOrder::RowMajor>;

namespace {

struct StrandCache {
	size_t pointOffset;
	size_t pointCount;
	std::vector<Float> prefixLength;
	Float length;

	StrandCache()
		: pointOffset(0), pointCount(0), length(0) {
	}
};

inline int clampIndex(Float value, int resolution) {
	return math::clamp(math::floorToInt(value), 0, resolution - 1);
}

inline Float clamp01(Float x) {
	return std::max((Float)0.0f, std::min((Float)1.0f, x));
}

bool evaluateStrandAtU(
	const std::vector<Point> &vertices,
	const StrandCache &strand,
	Float u,
	Point &p,
	Vector &tangent) {
	if (strand.pointCount < 2 || strand.length <= 1e-8f)
		return false;

	Float distance = clamp01(u) * strand.length;
	size_t segment = 0;
	while (segment + 1 < strand.pointCount - 1) {
		Float nextLength = strand.prefixLength[segment + 1] - strand.prefixLength[segment];
		if (distance <= strand.prefixLength[segment + 1] && nextLength > 1e-8f)
			break;
		segment++;
	}

	Float segStart = strand.prefixLength[segment];
	Float segEnd = strand.prefixLength[segment + 1];
	Float segLength = segEnd - segStart;
	size_t i0 = strand.pointOffset + segment;
	size_t i1 = i0 + 1;
	Point p0 = vertices[i0];
	Point p1 = vertices[i1];
	Vector dp = p1 - p0;
	if (segLength <= 1e-8f || dp.length() <= 1e-8f)
		return false;

	Float t = clamp01((distance - segStart) / segLength);
	p = p0 + dp * t;
	tangent = normalize(dp);
	return true;
}

Vector mulSymmetric3x3(Float m00, Float m01, Float m02, Float m11, Float m12, Float m22, const Vector &v) {
	return Vector(
		m00 * v.x + m01 * v.y + m02 * v.z,
		m01 * v.x + m11 * v.y + m12 * v.z,
		m02 * v.x + m12 * v.y + m22 * v.z);
}

Vector largestEigenvectorSymmetric3x3(Float m00, Float m01, Float m02, Float m11, Float m12, Float m22) {
	Vector seeds[4] = {
		Vector(1.0f, 0.0f, 0.0f),
		Vector(0.0f, 1.0f, 0.0f),
		Vector(0.0f, 0.0f, 1.0f),
		normalize(Vector(1.0f, 1.0f, 1.0f)),
	};

	Float bestValue = -std::numeric_limits<Float>::max();
	Vector bestVector(0.0f);
	for (const Vector &seed : seeds) {
		Vector v = seed;
		for (int i = 0; i < 16; ++i) {
			Vector mv = mulSymmetric3x3(m00, m01, m02, m11, m12, m22, v);
			Float len = mv.length();
			if (len <= 1e-8f)
				break;
			v = mv / len;
		}

		Float value = dot(v, mulSymmetric3x3(m00, m01, m02, m11, m12, m22, v));
		if (value > bestValue) {
			bestValue = value;
			bestVector = v;
		}
	}

	return bestVector.lengthSquared() > 1e-8f ? normalize(bestVector) : Vector(0.0f);
}

} // namespace

HairVoxel::HairVoxel(const Properties &props)
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

	dMultiplier = props.getFloat("dMultiplier", 1.5f);
	sampleDeltaAlongHair = props.getFloat("sampleDeltaAlongHair", 0.02f);
	if (dMultiplier <= 0.0f) {
		Log(EWarn, "HairVoxel: dMultiplier must be > 0. Resetting to 1.5.");
		dMultiplier = 1.5f;
	}
	if (sampleDeltaAlongHair <= 0.0f) {
		Log(EWarn, "HairVoxel: sampleDeltaAlongHair must be > 0. Resetting to 0.02.");
		sampleDeltaAlongHair = 0.02f;
	}

	outputSuffix = props.getString("outputSuffix", "");
	voxelData = HairVoxelGrid(voxelCountX, voxelCountY, voxelCountZ);
}

HairVoxel::HairVoxel(Stream *stream, InstanceManager *manager)
	: Precomputer(stream, manager) {
}

void HairVoxel::runPreBuild(Scene *scene) {
	Log(EInfo, "Precomputing HairVoxel.");
	ref<Timer> timer = new Timer();

	ref<BSplineHairShape> hairShape;
	auto shapes = scene->getShapes();
	for (auto &shape : shapes) {
		BSplineHairShape *candidate = dynamic_cast<BSplineHairShape *>(shape.get());
		if (!candidate)
			continue;
		if (hairShape) {
			Log(EError, "HairVoxel supports exactly one bSpline hair shape.");
			return;
		}
		hairShape = ref<BSplineHairShape>(candidate);
	}

	if (!hairShape) {
		Log(EError, "HairVoxel requires a bSpline hair shape in the scene.");
		return;
	}

	buildVoxelData(hairShape);

	std::string hairFileName = hairShape->getFileName();
	size_t dotPos = hairFileName.find_last_of('.');
	std::string hairName = (dotPos != std::string::npos) ? hairFileName.substr(0, dotPos) : hairFileName;
	serializeHairVoxelGrid(hairName, scene->getSourceFile());

	Log(EInfo, "Precomputing HairVoxel done. Time taken: %i ms.", timer->getMilliseconds());
}

Vector3i HairVoxel::getVoxelIndex(const Point &p) const {
	Vector3 rel = CwiseDiv((p - gridAABB.min), voxelSize);
	return Vector3i(
		math::clamp(math::floorToInt(rel.x), 0, voxelCountX - 1),
		math::clamp(math::floorToInt(rel.y), 0, voxelCountY - 1),
		math::clamp(math::floorToInt(rel.z), 0, voxelCountZ - 1));
}

Point HairVoxel::getVoxelCenter(int x, int y, int z) const {
	return gridAABB.min + CwiseMul(voxelSize, Vector((Float)x + 0.5f, (Float)y + 0.5f, (Float)z + 0.5f));
}

void HairVoxel::buildVoxelData(ref<BSplineHairShape> bSplineHairShape) {
	const std::vector<Point> &vertices = bSplineHairShape->getVertices();
	const std::vector<bool> &startsFiber = bSplineHairShape->getStartFiber();
	if (vertices.size() < 2 || startsFiber.size() < vertices.size()) {
		Log(EError, "HairVoxel: invalid strand vertex layout.");
		return;
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

	Float avgVoxelSize = (voxelSize.x + voxelSize.y + voxelSize.z) / 3.0f;
	Float d = dMultiplier * avgVoxelSize;
	Float fiberRadius = bSplineHairShape->getRadius();
	if (d <= 0.0f) {
		Log(EError, "HairVoxel: invalid gather radius.");
		return;
	}

	std::vector<StrandCache> strands;
	strands.reserve(32768);
	size_t pointOffset = 0;
	while (pointOffset < vertices.size()) {
		size_t nextOffset = pointOffset + 1;
		while (nextOffset < vertices.size() && !startsFiber[nextOffset])
			nextOffset++;

		size_t pointCount = nextOffset - pointOffset;
		if (pointCount >= 2) {
			StrandCache cache;
			cache.pointOffset = pointOffset;
			cache.pointCount = pointCount;
			cache.prefixLength.resize(pointCount, 0.0f);
			for (size_t i = 1; i < pointCount; ++i)
				cache.prefixLength[i] = cache.prefixLength[i - 1] + (vertices[pointOffset + i] - vertices[pointOffset + i - 1]).length();
			cache.length = cache.prefixLength.back();
			if (cache.length > 1e-8f)
				strands.push_back(cache);
		}

		pointOffset = nextOffset;
	}

	Log(EInfo, "HairVoxel: sampling %zu strands with sampleDeltaAlongHair %.5f.", strands.size(), sampleDeltaAlongHair);

	HairSampleVoxelGrid sampleGrid(voxelCountX, voxelCountY, voxelCountZ);
	size_t sampleCount = 0;
	for (uint32_t strandID = 0; strandID < (uint32_t)strands.size(); ++strandID) {
		const StrandCache &strand = strands[strandID];
		for (Float u = 0.0f; u < 1.0f; u += sampleDeltaAlongHair) {
			Point p;
			Vector tangent;
			if (!evaluateStrandAtU(vertices, strand, u, p, tangent))
				continue;

			Vector3 rel = CwiseDiv((p - gridAABB.min), voxelSize);
			int x = clampIndex(rel.x, voxelCountX);
			int y = clampIndex(rel.y, voxelCountY);
			int z = clampIndex(rel.z, voxelCountZ);
			sampleGrid(x, y, z).push_back({p, tangent, strandID});
			sampleCount++;
		}
	}

	voxelData = HairVoxelGrid(voxelCountX, voxelCountY, voxelCountZ);
	size_t occupiedVoxelCount = 0;
	Float maxDensity = 0.0f;
	std::vector<Vector> gatheredTangents;
	gatheredTangents.reserve(128);
	std::set<uint32_t> contributingStrands;

	for (int z = 0; z < voxelCountZ; ++z) {
		for (int y = 0; y < voxelCountY; ++y) {
			for (int x = 0; x < voxelCountX; ++x) {
				if (sampleGrid(x, y, z).empty())
					continue;

				Point voxelCenter = getVoxelCenter(x, y, z);
				Point searchMin = voxelCenter - Vector(d, d, d);
				Point searchMax = voxelCenter + Vector(d, d, d);
				Vector3 relMin = CwiseDiv((searchMin - gridAABB.min), voxelSize);
				Vector3 relMax = CwiseDiv((searchMax - gridAABB.min), voxelSize);
				int minX = clampIndex(relMin.x, voxelCountX);
				int minY = clampIndex(relMin.y, voxelCountY);
				int minZ = clampIndex(relMin.z, voxelCountZ);
				int maxX = clampIndex(relMax.x, voxelCountX);
				int maxY = clampIndex(relMax.y, voxelCountY);
				int maxZ = clampIndex(relMax.z, voxelCountZ);

				gatheredTangents.clear();
				contributingStrands.clear();
				Float c00 = 0.0f, c01 = 0.0f, c02 = 0.0f, c11 = 0.0f, c12 = 0.0f, c22 = 0.0f;

				for (int nz = minZ; nz <= maxZ; ++nz) {
					for (int ny = minY; ny <= maxY; ++ny) {
						for (int nx = minX; nx <= maxX; ++nx) {
							const std::vector<HairSamplePoint> &samples = sampleGrid(nx, ny, nz);
							for (const auto &sample : samples) {
								Float dist = (voxelCenter - sample.p).length();
								if (dist > d)
									continue;
								Float weight = 1.0f - dist / d;
								if (weight <= 1e-8f)
									continue;

								Vector t = normalize(sample.tangent);
								gatheredTangents.push_back(t);
								contributingStrands.insert(sample.strandID);

								c00 += weight * t.x * t.x;
								c01 += weight * t.x * t.y;
								c02 += weight * t.x * t.z;
								c11 += weight * t.y * t.y;
								c12 += weight * t.y * t.z;
								c22 += weight * t.z * t.z;
							}
						}
					}
				}

				if (gatheredTangents.empty())
					continue;

				Vector axis = largestEigenvectorSymmetric3x3(c00, c01, c02, c11, c12, c22);
				if (axis.lengthSquared() <= 1e-8f)
					continue;

				Float sumCos = 0.0f;
				Float sumCos2 = 0.0f;
				int n = 0;
				for (const Vector &t : gatheredTangents) {
					Float cosTheta = std::abs(dot(t, axis));
					sumCos += cosTheta;
					sumCos2 += cosTheta * cosTheta;
					++n;
				}

				Float variance = 0.0f;
				if (n >= 2)
					variance = std::max((Float)0.0f, (sumCos2 - sumCos * sumCos / (Float)n) / (Float)(n - 1));

				Float density = 2.0f * fiberRadius * (Float)contributingStrands.size() / ((Float)M_PI * d * d);
				HairVoxelT &cell = voxelData(x, y, z);
				cell.averageDirection = axis;
				cell.variance = variance;
				cell.perpAttenuation = density;
				maxDensity = std::max(maxDensity, density);
				occupiedVoxelCount++;
			}
		}
	}

	Log(EInfo, "HairVoxel: samples=%zu, occupiedVoxels=%zu, maxDensity=%f, dMultiplier=%f, sampleDelta=%f.",
		sampleCount, occupiedVoxelCount, maxDensity, dMultiplier, sampleDeltaAlongHair);
}

QuantizedDirection HairVoxel::quantizeDirection(const Vector3 &dir) const {
	if (dir.lengthSquared() < 1e-8f)
		return {255, 0};

	Vector3 nDir = normalize(dir);
	float theta = std::acos(std::max(-1.0f, std::min(1.0f, (float)nDir.z)));
	float phi = std::atan2((float)nDir.y, (float)nDir.x);
	if (phi < 0.0f)
		phi += 2.0f * (float)M_PI;

	uint8_t thetaByte = static_cast<uint8_t>(std::round((theta / (float)M_PI) * 254.0f));
	uint8_t phiByte = static_cast<uint8_t>(std::round((phi / (2.0f * (float)M_PI)) * 254.0f));
	return {thetaByte, phiByte};
}

void HairVoxel::writeVolHeader(std::ofstream &stream, int32_t encodingType, int32_t channels,
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

void HairVoxel::serializeHairVoxelGrid(const std::string &hairName, const fs::path &sceneFilePath) const {
	const Vector3i &res = voxelData.getResolution();
	if (voxelData.getSize() == 0) {
		Log(EError, "HairVoxel: Attempting to serialize an empty VoxelGrid.");
		return;
	}

	boost::filesystem::path outputDirectory = sceneFilePath.parent_path();
	std::stringstream ss;
	ss << hairName << "_" << res.x << "_" << res.y << "_" << res.z << outputSuffix;
	std::string baseFilename = ss.str();

	boost::filesystem::path dirPath = outputDirectory / (baseFilename + "_dir.vol");
	std::ofstream dirStream(dirPath.string(), std::ios::binary);
	if (!dirStream.is_open()) {
		Log(EError, "HairVoxel: Could not open file for writing: \"%s\"", dirPath.string().c_str());
		return;
	}
	writeVolHeader(dirStream, 4, 3, res, gridAABB);
	for (int z = 0; z < res.z; ++z)
		for (int y = 0; y < res.y; ++y)
			for (int x = 0; x < res.x; ++x) {
				QuantizedDirection qdir = quantizeDirection(voxelData(x, y, z).averageDirection);
				dirStream.write(reinterpret_cast<const char *>(&qdir.theta), 1);
				dirStream.write(reinterpret_cast<const char *>(&qdir.phi), 1);
			}
	dirStream.close();

	boost::filesystem::path varPath = outputDirectory / (baseFilename + "_var.vol");
	std::ofstream varStream(varPath.string(), std::ios::binary);
	if (!varStream.is_open()) {
		Log(EError, "HairVoxel: Could not open file for writing: \"%s\"", varPath.string().c_str());
		return;
	}
	writeVolHeader(varStream, 1, 1, res, gridAABB);
	for (int z = 0; z < res.z; ++z)
		for (int y = 0; y < res.y; ++y)
			for (int x = 0; x < res.x; ++x) {
				float variance = (float)voxelData(x, y, z).variance;
				varStream.write(reinterpret_cast<const char *>(&variance), sizeof(float));
			}
	varStream.close();

	boost::filesystem::path rawHairVolPath = outputDirectory / (baseFilename + "_atten_raw.hairvol");
	std::ofstream rawHairVolStream(rawHairVolPath.string(), std::ios::binary);
	if (!rawHairVolStream.is_open()) {
		Log(EError, "HairVoxel: Could not open file for writing: \"%s\"", rawHairVolPath.string().c_str());
		return;
	}
	writeVolHeader(rawHairVolStream, 1, 1, res, gridAABB);
	float maxDensity = 0.0f;
	for (int z = 0; z < res.z; ++z)
		for (int y = 0; y < res.y; ++y)
			for (int x = 0; x < res.x; ++x) {
				float density = (float)voxelData(x, y, z).perpAttenuation;
				rawHairVolStream.write(reinterpret_cast<const char *>(&density), sizeof(float));
				maxDensity = std::max(maxDensity, density);
			}
	rawHairVolStream.close();

	if (maxDensity <= 1e-9f) {
		maxDensity = 1.0f;
		Log(EWarn, "HairVoxel: Max density is close to 0, setting scale to 1.0.");
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

	boost::filesystem::path attenPath = outputDirectory / (baseFilename + "_atten_" + filenameSuffix + ".vol");
	std::ofstream attenStream(attenPath.string(), std::ios::binary);
	if (!attenStream.is_open()) {
		Log(EError, "HairVoxel: Could not open file for writing: \"%s\"", attenPath.string().c_str());
		return;
	}
	writeVolHeader(attenStream, 1, 1, res, gridAABB);
	for (int z = 0; z < res.z; ++z)
		for (int y = 0; y < res.y; ++y)
			for (int x = 0; x < res.x; ++x) {
				float normalizedDensity = (float)voxelData(x, y, z).perpAttenuation / maxDensity;
				normalizedDensity = std::max(0.0f, std::min(1.0f, normalizedDensity));
				attenStream.write(reinterpret_cast<const char *>(&normalizedDensity), sizeof(float));
			}
	attenStream.close();

	boost::filesystem::path xmlPath = outputDirectory / (hairName + outputSuffix + "_densityScale.xml");
	std::ofstream xmlStream(xmlPath.string());
	if (!xmlStream.is_open()) {
		Log(EError, "HairVoxel: Could not open XML file for writing: \"%s\"", xmlPath.string().c_str());
		return;
	}
	xmlStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
	xmlStream << "<scene version=\"0.6.0\">" << std::endl;
	xmlStream << "    <default name=\"DensityScale\" value=\"" << filenameSuffix << "\"/>" << std::endl;
	xmlStream << "    <default name=\"DensityValue\" value=\"" << rawFloatStr << "\"/>" << std::endl;
	xmlStream << "</scene>" << std::endl;
	xmlStream.close();

	Log(EInfo, "HairVoxel: Serialized voxel volumes to \"%s\" (max density %f).",
		outputDirectory.string().c_str(), maxDensity);
}

MTS_IMPLEMENT_CLASS_API(HairVoxel, false, Precomputer, HairVoxelAPI)
MTS_EXPORT_PLUGIN(HairVoxel, "Hair voxel representation");

MTS_NAMESPACE_END
