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

#if !defined(__HAIR_H)
#define __HAIR_H

#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

class HairKDTree;

/**
 * \brief Intersection shape structure for cylindrical hair
 * segments with miter joints. This class expects an ASCII file containing
 * a list of hairs made from segments. Each line should contain an X,
 * Y and Z coordinate separated by a space. An empty line indicates
 * the start of a new hair.
 */
class MTS_EXPORT HairShape : public Shape {
  public:
	/// Construct a new HairShape instance given a properties object
	HairShape(const Properties &props);

	/// Unserialize from a binary data stream
	HairShape(Stream *stream, InstanceManager *manager);

	/// Serialize to a binary data stream
	void serialize(Stream *stream, InstanceManager *manager) const;

	// =============================================================
	//! @{ \name Access the internal vertex data
	// =============================================================

	/// Return the list of vertices underlying the hair shape
	const std::vector<Point> &getVertices() const;

	/**
	 * Return a boolean list specifying whether a vertex
	 * marks the beginning of a new fiber
	 */
	const std::vector<bool> &getStartFiber() const;

	//! @}
	// =============================================================

	// =============================================================
	//! @{ \name Implementation of the \ref Shape interface
	// =============================================================

	bool rayIntersect(const Ray &ray, Float mint,
					  Float maxt, Float &t, void *temp) const;

	bool rayIntersect(const Ray &ray, Float mint, Float maxt) const;

	void fillIntersectionRecord(const Ray &ray,
								const void *temp, Intersection &its) const;

	ref<TriMesh> createTriMesh();

	const KDTreeBase<AABB> *getKDTree() const;

	AABB getAABB() const;

	Float getSurfaceArea() const;

	size_t getPrimitiveCount() const;

	size_t getEffectivePrimitiveCount() const;

	/**
	 * @brief Computes the next ray origin for a ray scattering off the hair surface.
	 *
	 * If the outgoing direction 'wo' points into the hair strand, this method
	 * advances the origin 'p' to the other side of the strand to simulate
	 * transmission and avoid self-intersection. Otherwise, it returns the
	 * original point 'p'.
	 *
	 * @param its The current intersection structure.
	 * @param wo The outgoing ray direction in world coordinates.
	 * @return The potentially offsetted ray origin for the next path segment.
	 */
	Point getNextRayOrigin(const Intersection &its, const Vector &wo) const;

	//! @}
	// =============================================================

	/// Return a human-readable representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
  private:
	void loadAsciiHair(const fs::path &path, const Transform &objectToWorld, Float dpThresh, Float reduction, ref<Random> random, std::vector<Point> &vertices, std::vector<bool> &vertexStartsFiber, size_t &nDegenerate, size_t &nSkipped);

	void loadBinaryHair(ref<FileStream> binaryStream, const Transform &objectToWorld, Float dpThresh, Float reduction, ref<Random> random, std::vector<Point> &vertices, std::vector<bool> &vertexStartsFiber, size_t &nDegenerate, size_t &nSkipped);

	void loadCemYuskelBinaryHair(const fs::path &path, const Transform &objectToWorld, Float dpThresh, Float reduction, ref<Random> random, std::vector<Point> &vertices, std::vector<bool> &vertexStartsFiber, size_t &nDegenerate, size_t &nSkipped);

	ref<HairKDTree> m_kdtree;
};

MTS_NAMESPACE_END

#endif /* __HAIR_H */
