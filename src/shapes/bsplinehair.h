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
#pragma once
#if !defined(__BSPLINEHAIR_H)
#define __BSPLINEHAIR_H

#include <mitsuba/render/shape.h>
#include <mitsuba/core/pmf.h>
MTS_NAMESPACE_BEGIN

#ifdef _WIN32
#ifdef BSPLINEHAIR_EXPORTS
#define BSPLINEHAIR_API __declspec(dllexport)
#else
#define BSPLINEHAIR_API __declspec(dllimport)
#endif
#else
#define BSPLINEHAIR_API __attribute__((visibility("default")))
#endif



// A structure to hold geometric information about a point on a hair curve.
struct CurvePointInfo {
	Point p;		// Position on the curve's centerline
	Vector tangent; // Normalized tangent vector (direction of the curve)
	Float radius;	// The radius of the hair at this point
};

/**
 * \brief Intersection shape structure for cylindrical hair
 * segments with miter joints. This class expects an ASCII file containing
 * a list of hairs made from segments. Each line should contain an X,
 * Y and Z coordinate separated by a space. An empty line indicates
 * the start of a new hair.
 */
class BSPLINEHAIR_API BSplineHairShape : public Shape {
  public:
	/// Construct a new HairShape instance given a properties object
	BSplineHairShape(const Properties &props);

	~BSplineHairShape();

	/// Unserialize from a binary data stream
	BSplineHairShape(Stream *stream, InstanceManager *manager);

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

	AABB getAABB() const;

	size_t getPrimitiveCount() const;

	size_t getEffectivePrimitiveCount() const;

	const std::vector<uint32_t> &getPrimitiveIndices() const;
	Float getRadius() const {
		return m_radius;
	}

	/**
	 * @brief Evaluates the geometric properties of the curve at a specific
	 * parametric location.
	 *
	 * This function provides information about the curve's centerline, including
	 * its position, Frenet frame (tangent, normal, bitangent), and radius.
	 *
	 * @param primID The ID of the B-spline segment to evaluate.
	 * @param u The parametric coordinate along the curve segment (from 0 to 1).
	 * @param info A reference to a CurvePointInfo struct that will be filled
	 *             with the geometric data.
	 * @return True if the evaluation was successful (e.g., primID is valid),
	 *         false otherwise.
	 */
	bool evaluatePrimitive(uint32_t primID, float u, CurvePointInfo &info) const;

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

	inline std::string getFileName() {
		return m_filename;
	}

	inline RTCGeometry getEmbreeGeometry() const override {
		return m_embree_geom;
	}
	//! @}
	// =============================================================

	/// Return a human-readable representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
  private:
	void loadAsciiHair(const fs::path &path, const Transform &objectToWorld, Float dpThresh, Float reduction, ref<Random> random, std::vector<Point> &vertices, std::vector<bool> &vertexStartsFiber, size_t &nDegenerate, size_t &nSkipped);

	void loadBinaryHair(ref<FileStream> binaryStream, const Transform &objectToWorld, Float dpThresh, Float reduction, ref<Random> random, std::vector<Point> &vertices, std::vector<bool> &vertexStartsFiber, size_t &nDegenerate, size_t &nSkipped);

	void loadCemYuskelBinaryHair(const fs::path &path, const Transform &objectToWorld, Float dpThresh, Float reduction, ref<Random> random, std::vector<Point> &vertices, std::vector<bool> &vertexStartsFiber, size_t &nDegenerate, size_t &nSkipped);

	void prepareDataForEmbreee(std::vector<Point> vertices, std::vector<bool> vertexStartsFiber);

	void buildEmbreeGeometry(RTCDevice device);

	// --- Internal Data Representation (Cleaned & World Space) ---

	// Stores all control points for all hair strands in a single flat list (World Space coordinates).
	std::vector<Point> m_original_vertices;

	// Parallel array acting as flags. 'true' means the corresponding vertex in m_original_vertices
	// is the first point (root) of a new hair strand.
	std::vector<bool> m_original_vertex_starts_fiber;

	// Member variables to store Embree objects
	RTCDevice m_embree_device = nullptr;
	RTCGeometry m_embree_geom = nullptr;
	RTCScene m_embree_scene = nullptr; // Each shape can have its own scene for intersection
	// --- Embree Hardware Acceleration Data (Formatted for Ray Tracing) ---

	// Flat buffer containing vertex data interleaved as (x, y, z, radius).
	// Size = total_vertices * 4.
	std::vector<float> m_embree_vertices;

	// Indices pointing to the start of each B-spline segment in m_embree_vertices.
	// Since cubic B-splines require 4 control points, index 'i' implies usage of points [i, i+1, i+2, i+3].
	std::vector<uint32_t> m_embree_indices;
	std::string m_filename;
};

MTS_NAMESPACE_END

#endif /* __BSPLINEHAIR_H */
