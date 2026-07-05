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

#include <mitsuba/render/volume.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/mmap.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

MTS_NAMESPACE_BEGIN

/*!\plugin{hairvolvolume}{Hair density volume data source}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Specifies the filename of the volume data file to be loaded
 *     }
 *     \parameter{sendData}{\Boolean}{
 *       When this parameter is set to \code{true}, the implementation will
 *       send all volume data to other network render nodes. Otherwise, they
 *       are expected to have access to an identical volume data file that can be
 *       mapped into memory. \default{\code{false}}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *         Optional linear transformation that should be applied to the data
 *     }
 *     \parameter{min, max}{\Point}{
 *         Optional parameter that can be used to re-scale the data so that
 *         it lies in the bounding box between \code{min} and \code{max}.
 *     }
 * }
 *
 * Reads the same VOL v3 binary layout as \pluginref{gridvolume}, but this plugin
 * is specialized for scalar density fields and returns the true maximum density
 * value stored in the file for Woodcock majorant construction.
 */
class HairVolDataSource : public VolumeDataSource {
public:
	enum EVolumeType {
		EFloat32 = 1,
		EFloat16 = 2,
		EUInt8 = 3,
		EQuantizedDirections = 4
	};

	HairVolDataSource(const Properties &props)
		: VolumeDataSource(props), m_data(nullptr), m_sendData(false),
		  m_volumeType(EFloat32), m_channels(1), m_stepSize(std::numeric_limits<Float>::infinity()),
		  m_maxDensity(0.0f) {
		m_volumeToWorld = props.getTransform("toWorld", Transform());

		if (props.hasProperty("min") && props.hasProperty("max")) {
			m_dataAABB.min = props.getPoint("min");
			m_dataAABB.max = props.getPoint("max");
		}

		m_sendData = props.getBoolean("sendData", false);
		loadFromFile(props.getString("filename"));
		configure();
	}

	HairVolDataSource(Stream *stream, InstanceManager *manager)
		: VolumeDataSource(stream, manager), m_data(nullptr), m_sendData(false),
		  m_volumeType(EFloat32), m_channels(1), m_stepSize(std::numeric_limits<Float>::infinity()),
		  m_maxDensity(0.0f) {
		m_volumeToWorld = Transform(stream);
		m_dataAABB = AABB(stream);
		m_sendData = stream->readBool();
		if (m_sendData) {
			m_volumeType = (EVolumeType)stream->readInt();
			m_res = Vector3i(stream);
			m_channels = stream->readInt();
			m_filename = stream->readString();
			size_t volumeSize = getVolumeSize();
			m_data = new uint8_t[volumeSize];
			stream->read(m_data, volumeSize);
		} else {
			fs::path filename = stream->readString();
			loadFromFile(filename);
		}
		configure();
	}

	virtual ~HairVolDataSource() {
		if (!m_mmap && m_data != nullptr)
			delete[] m_data;
	}

	size_t getVolumeSize() const {
		size_t nEntries = (size_t)m_res.x * (size_t)m_res.y * (size_t)m_res.z;
		switch (m_volumeType) {
			case EFloat32:
				return 4 * nEntries * m_channels;
			case EUInt8:
				return 1 * nEntries * m_channels;
			case EFloat16:
			case EQuantizedDirections:
			default:
				Log(EError, "HairVolDataSource only supports float32 or uint8 scalar volumes.");
				return 0;
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);

		m_volumeToWorld.serialize(stream);
		m_dataAABB.serialize(stream);
		stream->writeBool(m_sendData);

		if (m_sendData) {
			stream->writeInt(m_volumeType);
			m_res.serialize(stream);
			stream->writeInt(m_channels);
			stream->writeString(m_filename.string());
			stream->write(m_data, getVolumeSize());
		} else {
			stream->writeString(m_filename.string());
		}
	}

	void configure() {
		Vector extents(m_dataAABB.getExtents());
		m_worldToVolume = m_volumeToWorld.inverse();
		m_worldToGrid = Transform::scale(Vector(
			(m_res[0] - 1) / extents[0],
			(m_res[1] - 1) / extents[1],
			(m_res[2] - 1) / extents[2])) *
			Transform::translate(-Vector(m_dataAABB.min)) * m_worldToVolume;

		m_stepSize = std::numeric_limits<Float>::infinity();
		for (int i = 0; i < 3; ++i)
			m_stepSize = std::min(m_stepSize, 0.5f * extents[i] / (Float)(m_res[i] - 1));

		m_aabb.reset();
		for (int i = 0; i < 8; ++i)
			m_aabb.expandBy(m_volumeToWorld(m_dataAABB.getCorner(i)));

		computeCellMajorants();
		computeMaximumDensity();
	}

	void loadFromFile(const fs::path &filename) {
		m_filename = filename;
		fs::path resolved = Thread::getThread()->getFileResolver()->resolve(filename);
		m_mmap = new MemoryMappedFile(resolved);
		ref<MemoryStream> stream = new MemoryStream(m_mmap->getData(), m_mmap->getSize());
		stream->setByteOrder(Stream::ELittleEndian);

		char header[3];
		stream->read(header, 3);
		if (header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
			Log(EError, "Encountered an invalid volume data file (incorrect header identifier)");

		uint8_t version;
		stream->read(&version, 1);
		if (version != 3)
			Log(EError, "Encountered an invalid volume data file (incorrect file version)");

		int type = stream->readInt();
		int xres = stream->readInt(),
			yres = stream->readInt(),
			zres = stream->readInt();
		m_res = Vector3i(xres, yres, zres);
		m_channels = stream->readInt();

		switch (type) {
			case EFloat32:
			case EUInt8:
				if (m_channels != 1)
					Log(EError, "HairVolDataSource expects scalar density data (1 channel), got %i", m_channels);
				break;
			case EFloat16:
				Log(EError, "HairVolDataSource does not support float16 volumes.");
				break;
			case EQuantizedDirections:
				Log(EError, "HairVolDataSource does not support quantized direction volumes.");
				break;
			default:
				Log(EError, "Encountered a volume data file of unknown type (type=%i)", type);
		}

		m_volumeType = (EVolumeType)type;

		if (!m_dataAABB.isValid()) {
			Float xmin = stream->readSingle(),
				  ymin = stream->readSingle(),
				  zmin = stream->readSingle();
			Float xmax = stream->readSingle(),
				  ymax = stream->readSingle(),
				  zmax = stream->readSingle();
			m_dataAABB = AABB(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
		}

		Log(EDebug, "Mapped \"%s\" into memory: %ix%ix%i (%i channels), %s, %s",
			resolved.filename().string().c_str(), m_res.x, m_res.y, m_res.z, m_channels,
			memString(m_mmap->getSize()).c_str(), m_dataAABB.toString().c_str());

		m_data = (uint8_t *)m_mmap->getData() + 48;
	}

	Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = math::floorToInt(p.x),
				  y1 = math::floorToInt(p.y),
				  z1 = math::floorToInt(p.z),
				  x2 = x1 + 1, y2 = y1 + 1, z2 = z1 + 1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
			y2 >= m_res.y || z2 >= m_res.z)
			return 0.0f;

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
					_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		if (m_volumeType == EFloat32) {
			const float *floatData = (float *)m_data;
			const Float
				d000 = floatData[(z1 * m_res.y + y1) * m_res.x + x1],
				d001 = floatData[(z1 * m_res.y + y1) * m_res.x + x2],
				d010 = floatData[(z1 * m_res.y + y2) * m_res.x + x1],
				d011 = floatData[(z1 * m_res.y + y2) * m_res.x + x2],
				d100 = floatData[(z2 * m_res.y + y1) * m_res.x + x1],
				d101 = floatData[(z2 * m_res.y + y1) * m_res.x + x2],
				d110 = floatData[(z2 * m_res.y + y2) * m_res.x + x1],
				d111 = floatData[(z2 * m_res.y + y2) * m_res.x + x2];

			return ((d000 * _fx + d001 * fx) * _fy +
					(d010 * _fx + d011 * fx) * fy) *
					   _fz +
				   ((d100 * _fx + d101 * fx) * _fy +
					(d110 * _fx + d111 * fx) * fy) *
					   fz;
		}

		const Float inv255 = 1.0f / 255.0f;
		const Float
			d000 = m_data[(z1 * m_res.y + y1) * m_res.x + x1] * inv255,
			d001 = m_data[(z1 * m_res.y + y1) * m_res.x + x2] * inv255,
			d010 = m_data[(z1 * m_res.y + y2) * m_res.x + x1] * inv255,
			d011 = m_data[(z1 * m_res.y + y2) * m_res.x + x2] * inv255,
			d100 = m_data[(z2 * m_res.y + y1) * m_res.x + x1] * inv255,
			d101 = m_data[(z2 * m_res.y + y1) * m_res.x + x2] * inv255,
			d110 = m_data[(z2 * m_res.y + y2) * m_res.x + x1] * inv255,
			d111 = m_data[(z2 * m_res.y + y2) * m_res.x + x2] * inv255;

		return ((d000 * _fx + d001 * fx) * _fy +
				(d010 * _fx + d011 * fx) * fy) *
				   _fz +
			   ((d100 * _fx + d101 * fx) * _fy +
				(d110 * _fx + d111 * fx) * fy) *
				   fz;
	}

	bool supportsFloatLookups() const {
		return true;
	}

	bool supportsFloatLocalMajorant() const {
		return !m_cellMajorants.empty();
	}

	bool lookupFloatLocalMajorantSegment(const Point &p, const Vector &d,
		Float &majorant, Float &tMax) const {
		majorant = 0.0f;
		tMax = 0.0f;

		if (m_cellMajorants.empty())
			return false;

		const Point g = m_worldToGrid.transformAffine(p);
		const Vector dg = m_worldToGrid(d);

		const auto cellIndex = [](Float coord, Float dirComp) -> int {
			Float rel = coord;
			Float floorRel = std::floor(rel);
			Float frac = rel - floorRel;
			if (dirComp < 0.0f && frac < 1e-6f && rel > 0.0f)
				rel -= 1e-6f;
			return math::floorToInt(rel);
		};

		int ix = cellIndex(g.x, dg.x);
		int iy = cellIndex(g.y, dg.y);
		int iz = cellIndex(g.z, dg.z);
		if (ix < 0 || iy < 0 || iz < 0 || ix >= m_cellRes.x || iy >= m_cellRes.y || iz >= m_cellRes.z)
			return false;

		majorant = m_cellMajorants[(size_t)(ix + m_cellRes.x * (iy + m_cellRes.y * iz))];

		Float tExit = std::numeric_limits<Float>::infinity();
		const Float eps = 1e-8f;

		if (std::abs(dg.x) > eps) {
			Float next = (dg.x > 0.0f) ? (Float)(ix + 1) : (Float)ix;
			Float tx = (next - g.x) / dg.x;
			if (tx > eps)
				tExit = std::min(tExit, tx);
		}
		if (std::abs(dg.y) > eps) {
			Float next = (dg.y > 0.0f) ? (Float)(iy + 1) : (Float)iy;
			Float ty = (next - g.y) / dg.y;
			if (ty > eps)
				tExit = std::min(tExit, ty);
		}
		if (std::abs(dg.z) > eps) {
			Float next = (dg.z > 0.0f) ? (Float)(iz + 1) : (Float)iz;
			Float tz = (next - g.z) / dg.z;
			if (tz > eps)
				tExit = std::min(tExit, tz);
		}

		if (!std::isfinite(tExit))
			return false;

		tMax = tExit;
		return tMax > 0.0f;
	}

	bool supportsSpectrumLookups() const {
		return false;
	}

	bool supportsVectorLookups() const {
		return false;
	}

	Float getStepSize() const {
		return m_stepSize;
	}

	Float getMaximumFloatValue() const {
		return m_maxDensity;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HairVolVolume[" << endl
			<< "  res = " << m_res.toString() << "," << endl
			<< "  channels = " << m_channels << "," << endl
			<< "  maxDensity = " << m_maxDensity << "," << endl
			<< "  aabb = " << m_dataAABB.toString() << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	void computeCellMajorants() {
		m_cellRes = Vector3i(std::max(0, m_res.x - 1), std::max(0, m_res.y - 1), std::max(0, m_res.z - 1));
		size_t nCells = (size_t)m_cellRes.x * (size_t)m_cellRes.y * (size_t)m_cellRes.z;
		m_cellMajorants.clear();
		if (nCells == 0)
			return;

		m_cellMajorants.resize(nCells, 0.0f);

		if (m_volumeType == EFloat32) {
			const float *floatData = (const float *)m_data;
			for (int z = 0; z < m_cellRes.z; ++z) {
				for (int y = 0; y < m_cellRes.y; ++y) {
					for (int x = 0; x < m_cellRes.x; ++x) {
						Float m = 0.0f;
						for (int dz = 0; dz <= 1; ++dz) {
							for (int dy = 0; dy <= 1; ++dy) {
								for (int dx = 0; dx <= 1; ++dx) {
									size_t idx = (size_t)((z + dz) * m_res.y + (y + dy)) * (size_t)m_res.x + (size_t)(x + dx);
									Float v = (Float)floatData[idx];
									if (std::isfinite(v))
										m = std::max(m, v);
								}
							}
						}
						m_cellMajorants[(size_t)(x + m_cellRes.x * (y + m_cellRes.y * z))] = m;
					}
				}
			}
		} else if (m_volumeType == EUInt8) {
			const Float inv255 = 1.0f / 255.0f;
			for (int z = 0; z < m_cellRes.z; ++z) {
				for (int y = 0; y < m_cellRes.y; ++y) {
					for (int x = 0; x < m_cellRes.x; ++x) {
						uint8_t maxByte = 0;
						for (int dz = 0; dz <= 1; ++dz) {
							for (int dy = 0; dy <= 1; ++dy) {
								for (int dx = 0; dx <= 1; ++dx) {
									size_t idx = (size_t)((z + dz) * m_res.y + (y + dy)) * (size_t)m_res.x + (size_t)(x + dx);
									maxByte = std::max(maxByte, m_data[idx]);
								}
							}
						}
						m_cellMajorants[(size_t)(x + m_cellRes.x * (y + m_cellRes.y * z))] = (Float)maxByte * inv255;
					}
				}
			}
		}
	}

	void computeMaximumDensity() {
		m_maxDensity = 0.0f;
		size_t nEntries = (size_t)m_res.x * (size_t)m_res.y * (size_t)m_res.z;

		if (m_volumeType == EFloat32) {
			const float *floatData = (const float *)m_data;
			for (size_t i = 0; i < nEntries; ++i) {
				Float v = (Float)floatData[i];
				if (std::isfinite(v) && v > m_maxDensity)
					m_maxDensity = v;
			}
		} else if (m_volumeType == EUInt8) {
			uint8_t maxByte = 0;
			for (size_t i = 0; i < nEntries; ++i)
				maxByte = std::max(maxByte, m_data[i]);
			m_maxDensity = (Float)maxByte / 255.0f;
		}
	}

protected:
	fs::path m_filename;
	uint8_t *m_data;
	bool m_sendData;
	EVolumeType m_volumeType;
	Vector3i m_res;
	int m_channels;
	Transform m_worldToGrid;
	Transform m_worldToVolume;
	Transform m_volumeToWorld;
	Float m_stepSize;
	AABB m_dataAABB;
	ref<MemoryMappedFile> m_mmap;
	Float m_maxDensity;
	Vector3i m_cellRes;
	std::vector<Float> m_cellMajorants;
};

MTS_IMPLEMENT_CLASS_S(HairVolDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(HairVolDataSource, "Hair density volume data source");
MTS_NAMESPACE_END
