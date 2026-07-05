#pragma once
#if !defined(__VOXELGRID_H_)
#define __VOXELGRID_H_

#include <cstdint>
#include <mitsuba/mitsuba.h>
#include <vector>

MTS_NAMESPACE_BEGIN

#if !defined(NDEBUG)
#define DCHECK(condition)                                           \
	do {                                                            \
		if (!(condition)) {                                         \
			SLog(EError, "Debug Assertion Failed: \"%s\" in %s:%i", \
				 #condition, __FILE__, __LINE__);                   \
			std::abort();                                           \
		}                                                           \
	} while (false)
#else
#define DCHECK(condition) ((void)0)
#endif

/**
 * @brief Memory storage order for a VoxelGrid.
 */
enum class EStorageOrder {
	RowMajor, ///< Standard row-major layout (z, then y, then x).
	ZOrder	  ///< Z-order curve (Morton code) layout for 3D locality.
};

namespace internal {

/**
 * @brief Expands a 10-bit integer into 30 bits by interleaving two zero bits
 *        after each original bit. Used by 3D Morton encoding.
 */
inline uint32_t part1by2(uint32_t n) {
	n = (n | (n << 16)) & 0x030000FF;
	n = (n | (n << 8)) & 0x0300F00F;
	n = (n | (n << 4)) & 0x030C30C3;
	n = (n | (n << 2)) & 0x09249249;
	return n;
}

/**
 * @brief Encodes a 3D coordinate into a Morton code.
 */
inline uint32_t mortonEncode(uint32_t x, uint32_t y, uint32_t z) {
	return part1by2(x) | (part1by2(y) << 1) | (part1by2(z) << 2);
}

/**
 * @brief Inverse of Morton bit expansion. Extracts every third bit from a
 *        30-bit integer and compacts the result back to 10 bits.
 */
inline uint32_t compact1by2(uint32_t n) {
	n &= 0x09249249; // Keep bits 0, 3, 6, 9, ...
	n = (n ^ (n >> 2)) & 0x030C30C3;
	n = (n ^ (n >> 4)) & 0x0300F00F;
	n = (n ^ (n >> 8)) & 0x030000FF;
	n = (n ^ (n >> 16)) & 0x000003FF;
	return n;
}

/**
 * @brief Decodes a Morton code into a 3D integer coordinate.
 */
inline Vector3i mortonDecode(uint32_t m) {
	return Vector3i(
		compact1by2(m),
		compact1by2(m >> 1),
		compact1by2(m >> 2));
}

} // namespace internal

template <EStorageOrder VOrder>
struct IndexCalculator;

template <>
struct IndexCalculator<EStorageOrder::RowMajor> {
	static inline size_t getIndex(const Vector3i &p, const Vector3i &res) {
		return (size_t)p.z * res.y * res.x + (size_t)p.y * res.x + (size_t)p.x;
	}

	static inline Vector3i getVectorIndex(size_t index, const Vector3i &res) {
		size_t layerSize = (size_t)res.x * res.y;
		int z = (int)(index / layerSize);
		size_t rem = index % layerSize;
		int y = (int)(rem / res.x);
		int x = (int)(rem % res.x);
		return Vector3i(x, y, z);
	}
};

template <>
struct IndexCalculator<EStorageOrder::ZOrder> {
	static inline size_t getIndex(const Vector3i &p, const Vector3i &res) {
		return internal::mortonEncode(p.x, p.y, p.z);
	}

	static inline Vector3i getVectorIndex(size_t index, const Vector3i &res) {
		return internal::mortonDecode((uint32_t)index);
	}
};

/**
 * @brief Lightweight 3D voxel grid container.
 *
 * The grid stores values in a flat array and maps 3D integer coordinates to
 * linear indices through the selected storage order.
 */
template <typename T, EStorageOrder VOrder = EStorageOrder::RowMajor>
struct VoxelGrid {
  public:
	VoxelGrid()
		: m_resolution(0) {
	}

	explicit VoxelGrid(const Vector3i &resolution)
		: m_resolution(resolution) {
		if (m_resolution.x < 0 || m_resolution.y < 0 || m_resolution.z < 0) {
			SLog(EError, "VoxelGrid resolution must be non-negative!");
			m_resolution = Vector3i(0);
		}
		size_t size = (size_t)m_resolution.x * m_resolution.y * m_resolution.z;
		m_data.resize(size);
	}

	VoxelGrid(const Vector3i &resolution, const T &defaultValue)
		: m_resolution(resolution) {
		if (m_resolution.x < 0 || m_resolution.y < 0 || m_resolution.z < 0) {
			SLog(EError, "VoxelGrid resolution must be non-negative!");
			m_resolution = Vector3i(0);
		}
		size_t size = (size_t)m_resolution.x * m_resolution.y * m_resolution.z;
		m_data.assign(size, defaultValue);
	}

	VoxelGrid(int resX, int resY, int resZ)
		: VoxelGrid(Vector3i(resX, resY, resZ)) {
	}

	T &operator[](size_t index) {
		DCHECK(index < m_data.size());
		return m_data[index];
	}

	const T &operator[](size_t index) const {
		DCHECK(index < m_data.size());
		return m_data[index];
	}

	T &at_linear(size_t index) {
		if (index >= m_data.size())
			throw std::out_of_range("VoxelGrid linear index out of range");
		return m_data[index];
	}

	using iterator = typename std::vector<T>::iterator;
	using const_iterator = typename std::vector<T>::const_iterator;

	iterator begin() {
		return m_data.begin();
	}

	iterator end() {
		return m_data.end();
	}

	const_iterator begin() const {
		return m_data.begin();
	}

	const_iterator end() const {
		return m_data.end();
	}

	const T &operator()(int x, int y, int z) const {
		DCHECK(x >= 0 && x < m_resolution.x &&
			   y >= 0 && y < m_resolution.y &&
			   z >= 0 && z < m_resolution.z);
		size_t index = IndexCalculator<VOrder>::getIndex(Vector3i(x, y, z), m_resolution);
		return m_data[index];
	}

	T &operator()(int x, int y, int z) {
		DCHECK(x >= 0 && x < m_resolution.x &&
			   y >= 0 && y < m_resolution.y &&
			   z >= 0 && z < m_resolution.z);
		size_t index = IndexCalculator<VOrder>::getIndex(Vector3i(x, y, z), m_resolution);
		return m_data[index];
	}

	const T &operator()(const Vector3i &p) const {
		DCHECK(p.x >= 0 && p.x < m_resolution.x &&
			   p.y >= 0 && p.y < m_resolution.y &&
			   p.z >= 0 && p.z < m_resolution.z);
		size_t index = IndexCalculator<VOrder>::getIndex(p, m_resolution);
		return m_data[index];
	}

	T &operator()(const Vector3i &p) {
		DCHECK(p.x >= 0 && p.x < m_resolution.x &&
			   p.y >= 0 && p.y < m_resolution.y &&
			   p.z >= 0 && p.z < m_resolution.z);
		size_t index = IndexCalculator<VOrder>::getIndex(p, m_resolution);
		return m_data[index];
	}

	const Vector3i &getResolution() const {
		return m_resolution;
	}

	size_t getSize() const {
		return (size_t)m_resolution.x * m_resolution.y * m_resolution.z;
	}

	void resize(const Vector3i &newResolution) {
		m_resolution = newResolution;
		if (m_resolution.x < 0 || m_resolution.y < 0 || m_resolution.z < 0) {
			m_resolution = Vector3i(0);
		}
		size_t size = (size_t)m_resolution.x * m_resolution.y * m_resolution.z;
		m_data.clear();
		m_data.resize(size);
	}

	const T *getData() const {
		return m_data.data();
	}

	T *getData() {
		return m_data.data();
	}

	void swap(VoxelGrid &other) noexcept {
		m_data.swap(other.m_data);
		std::swap(m_resolution, other.m_resolution);
	}

	inline size_t getLinearIndex(const Vector3i &p) const {
		DCHECK(p.x >= 0 && p.x < m_resolution.x &&
			   p.y >= 0 && p.y < m_resolution.y &&
			   p.z >= 0 && p.z < m_resolution.z);
		return IndexCalculator<VOrder>::getIndex(p, m_resolution);
	}

	inline Vector3i getVectorIndex(size_t index) const {
		DCHECK(index < getSize());
		return IndexCalculator<VOrder>::getVectorIndex(index, m_resolution);
	}

  private:
	std::vector<T> m_data;
	Vector3i m_resolution;
};

template <typename T, EStorageOrder VOrder>
void swap(VoxelGrid<T, VOrder> &a, VoxelGrid<T, VOrder> &b) noexcept {
	a.swap(b);
}

MTS_NAMESPACE_END

#endif
