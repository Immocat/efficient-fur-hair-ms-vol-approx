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

#include "bsplinehair.h"
#include <cyHair/cyHairFile.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/sahkdtree3.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/trimesh.h>

#define MTS_HAIR_USE_FANCY_CLIPPING 1

MTS_NAMESPACE_BEGIN

/*!\plugin{hair}{Hair intersection shape}
 * \order{11}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the hair data file that should be loaded
 *     }
 *     \parameter{radius}{\Float}{
 *       Radius of the hair segments in world-space units
 *       \default{0.025, which assumes that the scene
 *       is modeled in millimeters.}.
 *     }
 *     \parameter{angleThreshold}{\Float}{
 *       For performance reasons, the plugin will merge adjacent hair
 *       segments when the angle of their tangent directions is below
 *       than this value (in degrees). \default{1}.
 *     }
 *     \parameter{reduction}{\Float}{
 *       When the reduction ratio is set to a value between zero and one, the hair
 *       plugin stochastically culls this portion of the input data (where
 *       1 corresponds to removing all hairs). To approximately preserve the
 *       appearance in renderings, the hair radius is enlarged (see Cook et al.
 *       \cite{Cook2007Stochastic}). This parameter is convenient for fast
 *       previews. \default{0, i.e. all geometry is rendered}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *        Specifies an optional linear object-to-world transformation.
 *        Note that non-uniform scales are not permitted!
 *        \default{none, i.e. object space $=$ world space}
 *     }
 * }
 * \renderings{
 *     \centering
 *     \fbox{\includegraphics[width=6cm]{images/shape_hair}}\hspace{4.5cm}
 *     \caption{A close-up of the hair shape rendered with a diffuse
 *     scattering model (an actual hair scattering model will
 *     be needed for realistic appearance)}
 * }
 * The plugin implements a space-efficient acceleration structure for
 * hairs made from many straight cylindrical hair segments with miter
 * joints. The underlying idea is that intersections with straight cylindrical
 * hairs can be found quite efficiently, and curved hairs are easily
 * approximated using a series of such segments.
 *
 * The plugin supports two different input formats: a simple (but not
 * particularly efficient) ASCII format containing the coordinates of a
 * hair vertex on every line. An empty line marks the beginning of a
 * new hair. The following snippet is an example of this format:\newpage
 * \begin{xml}
 * .....
 * -18.5498 -21.7669 22.8138
 * -18.6358 -21.3581 22.9262
 * -18.7359 -20.9494 23.0256
 *
 * -30.6367 -21.8369 6.78397
 * -30.7289 -21.4145 6.76688
 * -30.8226 -20.9933 6.73948
 * .....
 * \end{xml}
 *
 * There is also a binary format, which starts with the identifier
 * ``\texttt{BINARY\_HAIR}'' (11 bytes), followed by the number of
 * vertices as a 32-bit little endian integer.
 * The remainder of the file consists of the vertex positions stored as
 * single-precision XYZ coordinates (again in little-endian byte ordering).
 * To mark the beginning of a new hair strand, a single $+\infty$ floating
 * point value can be inserted between the vertex data.
 */

BSplineHairShape::BSplineHairShape(const Properties &props)
	: Shape(props) {
	m_isCurve = true;

	fs::path path = Thread::getThread()->getFileResolver()->resolve(
		props.getString("filename"));

	m_filename = props.getString("filename");

	m_radius = props.getFloat("radius", 0.025f);
	/* Skip segments, whose tangent differs by less than one degree
	   compared to the previous one */
	Float angleThreshold = degToRad(props.getFloat("angleThreshold", 1.0f));
	Float dpThresh = std::cos(angleThreshold);

	/* When set to a value n>1, the hair shape object will reduce
	   the input by only loading every n-th hair */
	Float reduction = props.getFloat("reduction", 0);
	if (reduction < 0 || reduction >= 1) {
		Log(EError, "The 'reduction' parameter must have a value in [0, 1)!");
	} else if (reduction > 0) {
		Float correction = 1.0f / (1 - reduction);
		Log(EDebug, "Reducing the amount of geometry by %.2f%%, scaling radii by %f.",
			reduction * 100, correction);
		m_radius *= correction;
	}
	m_thickness = m_radius * 2;

	/* Object-space -> World-space transformation */
	Transform objectToWorld = props.getTransform("toWorld", Transform());
	// m_radius *= objectToWorld(Vector(0, 0, 1)).length();

	Log(EInfo, "Loading hair geometry from \"%s\" ..", path.filename().string().c_str());
	ref<Timer> timer = new Timer();

	std::vector<Point> vertices;
	std::vector<bool> vertexStartsFiber;
	size_t nDegenerate = 0, nSkipped = 0;
	ref<Random> random = new Random();

	ref<FileStream> fs = new FileStream(path, FileStream::EReadOnly);
	fs->setByteOrder(Stream::ELittleEndian);

	const char *binaryHeader = "BINARY_HAIR";
	const char *cemYukselBinaryHeader = "HAIR";
	char temp[11];
	enum class HairFormat { ASCII,
							BINARY_HAIR,
							CEM_YUKSEL_HAIR };
	HairFormat format = HairFormat::ASCII;

	fs->seek(0);
	if (fs->getSize() >= 4) {
		fs->read(temp, 4);
		if (memcmp(temp, cemYukselBinaryHeader, 4) == 0) {
			format = HairFormat::CEM_YUKSEL_HAIR;
		}
	}

	if (format == HairFormat::ASCII && fs->getSize() >= 11) {
		fs->seek(0);
		fs->read(temp, 11);
		if (memcmp(temp, binaryHeader, 11) == 0) {
			format = HairFormat::BINARY_HAIR;
		}
	}
	fs->seek(0);

	switch (format) {
	case HairFormat::CEM_YUKSEL_HAIR:
		Log(EInfo, "Detected 'HAIR' binary format.");
		loadCemYuskelBinaryHair(path, objectToWorld, dpThresh, reduction, random,
								vertices, vertexStartsFiber, nDegenerate, nSkipped);
		break;

	case HairFormat::BINARY_HAIR:
		Log(EInfo, "Detected 'BINARY_HAIR' format.");
		loadBinaryHair(fs, objectToWorld, dpThresh, reduction, random,
					   vertices, vertexStartsFiber, nDegenerate, nSkipped);
		break;

	case HairFormat::ASCII:
	default:
		Log(EInfo, "Assuming ASCII format.");
		fs = nullptr;
		fs->close();
		loadAsciiHair(path, objectToWorld, dpThresh, reduction, random,
					  vertices, vertexStartsFiber, nDegenerate, nSkipped);
		break;
	}

	if (nDegenerate > 0)
		Log(EInfo, "Encountered " SIZE_T_FMT " degenerate segments!", nDegenerate);
	if (nSkipped > 0)
		Log(EInfo, "Skipped " SIZE_T_FMT " segments.", nSkipped);
	Log(EInfo, "Done (took %i ms)", timer->getMilliseconds());

	vertexStartsFiber.push_back(true);

	m_original_vertices = vertices;
	m_original_vertex_starts_fiber = vertexStartsFiber;

	prepareDataForEmbreee(vertices, vertexStartsFiber);
	m_embree_device = rtcNewDevice(NULL); // NULL creates a default device
	// Handle potential errors
	rtcSetDeviceErrorFunction(
		m_embree_device, [](void *userPtr, RTCError error, const char *str) { Log(EError, "Embree error: %s", str); }, nullptr);
	// Call our new method to build the Embree structures
	if (!m_embree_indices.empty()) {
		buildEmbreeGeometry(m_embree_device);
	}

	auto aabb = this->getAABB();
	std::cout << aabb.toString() << std::endl;
}

// ----- Helper function for loading the BINARY_HAIR format -----
void BSplineHairShape::loadBinaryHair(
	// Input parameters
	ref<FileStream> binaryStream,
	const Transform &objectToWorld,
	Float dpThresh,
	Float reduction,
	ref<Random> random,

	// Output parameters (passed by reference)
	std::vector<Point> &vertices,
	std::vector<bool> &vertexStartsFiber,
	size_t &nDegenerate,
	size_t &nSkipped) {
	// These variables hold the state during parsing
	Vector tangent(0.0f);
	Point p, lastP(0.0f);
	bool ignore = false;

	// Skip the header
	char header[11];
	binaryStream->read(header, 11);

	size_t vertexCount = binaryStream->readUInt();
	Log(EInfo, "Loading " SIZE_T_FMT " hair vertices (binary format)..", vertexCount);
	vertices.reserve(vertexCount);
	vertexStartsFiber.reserve(vertexCount);

	bool newFiber = true;
	size_t verticesRead = 0;

	while (verticesRead != vertexCount) {
		Float value = binaryStream->readSingle();
		// Check for the new-hair marker (+infinity)
		if (std::isinf(value)) {
			// The next 3 floats are the vertex coordinates
			p.x = binaryStream->readSingle();
			p.y = binaryStream->readSingle();
			p.z = binaryStream->readSingle();
			newFiber = true;
			if (reduction > 0)
				ignore = random->nextFloat() < reduction;
		} else {
			p.x = value;
			p.y = binaryStream->readSingle();
			p.z = binaryStream->readSingle();
		}

		p = objectToWorld(p);
		verticesRead++;

		if (ignore) {
			++nSkipped;
		} else if (newFiber) {
			vertices.push_back(p);
			vertexStartsFiber.push_back(true); // Mark as start of a new fiber
			lastP = p;
			tangent = Vector(0.0f);
		} else if (p != lastP) {
			if (tangent.isZero()) {
				vertices.push_back(p);
				vertexStartsFiber.push_back(false);
				tangent = normalize(p - lastP);
				lastP = p;
			} else {
				Vector nextTangent = normalize(p - lastP);
				if (dot(nextTangent, tangent) > dpThresh) {
					// Tangent difference is too small, overwrite previous vertex
					tangent = normalize(p - vertices[vertices.size() - 2]);
					vertices[vertices.size() - 1] = p;
					++nSkipped;
				} else {
					vertices.push_back(p);
					vertexStartsFiber.push_back(false);
					tangent = nextTangent;
				}
				lastP = p;
			}
		} else {
			nDegenerate++;
		}
		newFiber = false;
	}
}

// ----- Helper function for loading the ASCII hair format -----
void BSplineHairShape::loadAsciiHair(
	// Input parameters
	const fs::path &path,
	const Transform &objectToWorld,
	Float dpThresh,
	Float reduction,
	ref<Random> random,

	// Output parameters (passed by reference)
	std::vector<Point> &vertices,
	std::vector<bool> &vertexStartsFiber,
	size_t &nDegenerate,
	size_t &nSkipped) {
	// These variables hold the state during parsing
	Vector tangent(0.0f);
	Point p, lastP(0.0f);
	bool ignore = false;
	bool newFiber = true;
	std::string line;

	fs::ifstream is(path);
	if (is.fail())
		Log(EError, "Could not open \"%s\"!", path.string().c_str());

	Log(EInfo, "Loading hair vertices (ASCII format)..");

	while (std::getline(is, line)) {
		// Skip comments
		if (line.length() > 0 && line[0] == '#') {
			continue;
		}

		std::istringstream iss(line);
		iss >> p.x >> p.y >> p.z;

		if (!iss.fail()) {
			// Successfully parsed a vertex
			p = objectToWorld(p);
			if (ignore) {
				++nSkipped;
			} else if (newFiber) {
				vertices.push_back(p);
				vertexStartsFiber.push_back(true); // Mark as start of a new fiber
				lastP = p;
				tangent = Vector(0.0f);
			} else if (p != lastP) {
				if (tangent.isZero()) {
					vertices.push_back(p);
					vertexStartsFiber.push_back(false);
					tangent = normalize(p - lastP);
					lastP = p;
				} else {
					Vector nextTangent = normalize(p - lastP);
					if (dot(nextTangent, tangent) > dpThresh) {
						// Tangent difference is too small, overwrite previous vertex
						tangent = normalize(p - vertices[vertices.size() - 2]);
						vertices[vertices.size() - 1] = p;
						++nSkipped;
					} else {
						vertices.push_back(p);
						vertexStartsFiber.push_back(false);
						tangent = nextTangent;
					}
					lastP = p;
				}
			} else {
				nDegenerate++;
			}
			newFiber = false;
		} else {
			// An empty or non-parseable line marks the beginning of a new hair
			newFiber = true;
			if (reduction > 0)
				ignore = random->nextFloat() < reduction;
		}
	}
}

// ----- Helper function for loading the Cem Yuskel hair format -----
// https://www.cemyuksel.com/research/hairmodels/
void BSplineHairShape::loadCemYuskelBinaryHair(
	// Input parameters
	const fs::path &path,
	const Transform &objectToWorld,
	Float dpThresh,
	Float reduction,
	ref<Random> random,

	// Output parameters (passed by reference)
	std::vector<Point> &vertices,
	std::vector<bool> &vertexStartsFiber,
	size_t &nDegenerate,
	size_t &nSkipped) {

	// Helper object to load the .hair file
	cyHairFile hairFile;

	// Load the hair data from the specified file path
	int hairCount = hairFile.LoadFromFile(path.string().c_str());
	if (hairCount <= 0) {
		Log(EError, "Failed to load hair file: %s", path.string());
		return;
	}

	// Get the header which contains metadata about the hair file
	const cyHairFile::Header &header = hairFile.GetHeader();

	// Log the number of strands and vertices being loaded
	Log(EInfo, "Loading " SIZE_T_FMT " hair strands, " SIZE_T_FMT " vertices (Cem Yuksel .hair format)..",
		(size_t)header.hair_count, (size_t)header.point_count);

	// --- Check for unsupported per-point attributes and issue warnings ---
	if (header.arrays & _CY_HAIR_FILE_THICKNESS_BIT) {
		Log(EWarn, "The hair file contains per-point m_thickness data, which is not supported by Mitsuba. "
				   "The default m_thickness from the shape's material will be used instead.");
	}
	if (header.arrays & _CY_HAIR_FILE_TRANSPARENCY_BIT) {
		Log(EWarn, "The hair file contains per-point transparency data, which is not supported by Mitsuba. "
				   "The default transparency from the shape's material will be used instead.");
	}
	if (header.arrays & _CY_HAIR_FILE_COLORS_BIT) {
		Log(EWarn, "The hair file contains per-point color data, which is not supported by Mitsuba. "
				   "The default color from the shape's material will be used instead.");
	}

	// Get pointers to the main data arrays from the file
	const float *pointsArray = hairFile.GetPointsArray();
	const unsigned short *segmentsArray = hairFile.GetSegmentsArray();

	// Reserve memory in the output vectors to avoid multiple reallocations, improving performance.
	vertices.reserve(header.point_count);
	vertexStartsFiber.reserve(header.point_count);

	if (!(header.arrays & _CY_HAIR_FILE_POINTS_BIT) || pointsArray == nullptr) {
		Log(EError, "Points array missing in hair file, cannot load.");
		return;
	}

	// Initialize counters
	nDegenerate = 0;
	nSkipped = 0;
	// This index tracks the current position in the flat 'points_array'
	size_t currentPointIndex = 0;

	// Iterate over each hair strand
	for (uint32_t i = 0; i < header.hair_count; ++i) {
		// Determine the number of segments for the current strand.
		// If the segments_array is not present, use the default value from the header.
		unsigned int segmentCount = segmentsArray ? segmentsArray[i] : header.d_segments;
		// The number of points in a strand is always one more than the number of segments.
		unsigned int pointsInStrand = segmentCount + 1;

		// Randomly skip a fraction of hairs based on the 'reduction' parameter
		if (reduction > 0 && random->nextFloat() < reduction) {
			currentPointIndex += pointsInStrand;
			nSkipped += pointsInStrand;
			continue; // Skip to the next strand
		}

		// These variables hold the state for processing a single fiber
		Vector tangent(0.0f);
		Point p, lastP(0.0f);
		bool isNewFiber = true;

		// Iterate over each point within the current hair strand
		for (uint32_t j = 0; j < pointsInStrand; ++j) {
			// Read the (x, y, z) coordinates for the current point
			p.x = pointsArray[currentPointIndex * 3 + 0];
			p.y = pointsArray[currentPointIndex * 3 + 1];
			p.z = pointsArray[currentPointIndex * 3 + 2];

			// Transform the point from object space to world space
			p = objectToWorld(p);

			if (isNewFiber) {
				vertices.push_back(p);
				vertexStartsFiber.push_back(true); // Mark this vertex as the start of a new fiber
				lastP = p;
				tangent = Vector(0.0f); // Reset tangent for the new fiber
				isNewFiber = false;
			} else if (p != lastP) {
				if (tangent.isZero()) {
					// This is the second point of the fiber, initialize the first tangent
					vertices.push_back(p);
					vertexStartsFiber.push_back(false);
					tangent = normalize(p - lastP);
					lastP = p;
				} else {
					// For subsequent points, check if they are nearly collinear with the previous point
					Vector nextTangent = normalize(p - lastP);
					// if (dot(nextTangent, tangent) > dpThresh) {
					//	// The change in tangent is too small, indicating a nearly straight line.
					//	// Overwrite the previous vertex to reduce redundant points.
					//	tangent = normalize(p - vertices[vertices.size() - 2]);
					//	vertices.back() = p;
					//	++nSkipped;
					// } else {
					//	// The point is not collinear, add it to the fiber
					//	vertices.push_back(p);
					//	vertexStartsFiber.push_back(false);
					//	tangent = nextTangent;
					// }
					vertices.push_back(p);
					vertexStartsFiber.push_back(false);
					tangent = nextTangent;
					lastP = p;
				}
			} else {
				// The current point is identical to the previous one (degenerate segment)
				++nDegenerate;
			}

			// Move to the next point in the global points array
			currentPointIndex++;
		}
	}
}

void BSplineHairShape::prepareDataForEmbreee(std::vector<Point> vertices, std::vector<bool> vertexStartsFiber) {
	Log(EInfo, "Generating B-spline segments for Embree..");
	ref<Timer> timer = new Timer();

	// 1. Prepare Embree-compatible buffers
	// Vertex buffer stores (x, y, z, radius) for each control point
	m_embree_vertices.reserve(vertices.size() * 4);
	// Index buffer stores the starting index of each 4-point B-spline segment
	m_embree_indices.reserve(vertices.size()); // A rough over-estimation is fine

	uint32_t embree_vertex_offset = 0;

	// 2. Iterate through all loaded vertices to process one strand at a time
	for (size_t i = 0; i < vertices.size(); ++i) {
		// Find the start of the current strand
		size_t strand_start_idx = i;

		// Find the end of the current strand by looking for the next 'true' marker
		size_t strand_end_idx = strand_start_idx;
		while (strand_end_idx + 1 < vertices.size() && !vertexStartsFiber[strand_end_idx + 1]) {
			strand_end_idx++;
		}

		size_t num_points_in_strand = strand_end_idx - strand_start_idx + 1;

		// A. Populate the Embree vertex buffer for this single strand
		for (size_t j = strand_start_idx; j <= strand_end_idx; ++j) {
			const Point &p = vertices[j];
			m_embree_vertices.push_back(p.x);
			m_embree_vertices.push_back(p.y);
			m_embree_vertices.push_back(p.z);
			m_embree_vertices.push_back(m_radius); // Add the uniform radius
		}

		// B. Create B-spline segments using the sliding window approach from the reference
		// A renderable B-spline segment requires 4 control points.
		if (num_points_in_strand >= 4) {
			// A strand with N points can produce N-3 overlapping segments.
			// This is equivalent to the reference's `while (t < s - 2)` where s = num_segments.
			for (size_t t = 0; t <= num_points_in_strand - 4; ++t) {
				// The index buffer points to the *first* control point of the segment
				// within the global Embree vertex buffer.
				m_embree_indices.push_back(embree_vertex_offset + t);
			}
		}

		// C. Update the global vertex offset for the next strand
		embree_vertex_offset += num_points_in_strand;

		// D. Jump the main loop counter to the end of the strand we just processed
		i = strand_end_idx;
	}

	// Cleanup the extra marker we added
	vertexStartsFiber.pop_back();

	Log(EInfo, "Created " SIZE_T_FMT " B-spline segments from " SIZE_T_FMT " vertices (took %i ms)",
		m_embree_indices.size(), m_embree_vertices.size() / 4, timer->getMilliseconds());
}

void BSplineHairShape::buildEmbreeGeometry(RTCDevice device) {
	// 1. Create the Geometry Object
	// This tells Embree what kind of primitive we are creating.
	// We use RTC_GEOMETRY_TYPE_ROUND_BSPLINE_CURVE for B-splines with thickness.
	m_embree_geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ROUND_BSPLINE_CURVE);

	// 2. Set the Vertex Buffer
	// This is where we hand our (x, y, z, radius) data to Embree.
	// We use `rtcSetSharedGeometryBuffer` to avoid copying the data. Embree will
	// point directly to our vector's memory.
	rtcSetSharedGeometryBuffer(
		m_embree_geom,				 // The geometry object
		RTC_BUFFER_TYPE_VERTEX,		 // We are specifying vertices
		0,							 // The buffer slot (usually 0)
		RTC_FORMAT_FLOAT4,			 // The format is 4 floats (x,y,z,r)
		m_embree_vertices.data(),	 // Pointer to the raw data
		0,							 // Byte offset to the first element
		sizeof(float[4]),			 // Stride: bytes between consecutive vertices
		m_embree_vertices.size() / 4 // Total number of vertices
	);

	// 3. Set the Index Buffer
	// This buffer tells Embree how to form curve segments from the vertex buffer.
	// Each integer is the starting index of a 4-point control polygon.
	rtcSetSharedGeometryBuffer(
		m_embree_geom,			 // The geometry object
		RTC_BUFFER_TYPE_INDEX,	 // We are specifying indices
		0,						 // The buffer slot (usually 0)
		RTC_FORMAT_UINT,		 // The format is unsigned int
		m_embree_indices.data(), // Pointer to the raw data
		0,						 // Byte offset to the first element
		sizeof(uint32_t),		 // Stride: bytes between consecutive indices
		m_embree_indices.size()	 // Total number of indices (and thus segments)
	);

	// 4. Commit the Geometry
	// This tells Embree that we are done defining this piece of geometry.
	// Embree will now build its internal data structures for it.
	rtcCommitGeometry(m_embree_geom);

	// 5. Create a Scene to hold the Geometry
	// A scene is a container for one or more geometries. Rays are traced against scenes.
	m_embree_scene = rtcNewScene(device);

	// 6. Attach the Geometry to the Scene
	// This adds our hair geometry to the scene and gives it a unique ID (in this case, 0).
	// If you had multiple geometries, they would get IDs 1, 2, 3...
	// This ID is returned in the ray hit information.
	rtcAttachGeometry(m_embree_scene, m_embree_geom);

	// 7. Commit the Scene
	// This tells Embree to build the top-level acceleration structure (BVH) for the scene.
	// After this, the scene is ready for ray tracing.
	rtcCommitScene(m_embree_scene);
}

BSplineHairShape::~BSplineHairShape() {
	if (m_embree_geom) {
		// Geometries attached to a scene are automatically released
		// when the scene is released, but it's good practice to do it explicitly.
		// First detach, then release.
		rtcDetachGeometry(m_embree_scene, 0); // Detach geometry with ID 0
		rtcReleaseGeometry(m_embree_geom);
	}
	if (m_embree_scene) {
		rtcReleaseScene(m_embree_scene);
	}

	if (m_embree_device) {
		rtcReleaseDevice(m_embree_device);
	}
}

BSplineHairShape::BSplineHairShape(Stream *stream, InstanceManager *manager)
	: Shape(stream, manager) {
	Float radius = stream->readFloat();
	size_t vertexCount = stream->readSize();

	std::vector<Point> vertices(vertexCount);
	std::vector<bool> vertexStartsFiber(vertexCount + 1);
	stream->readFloatArray((Float *)&vertices[0], vertexCount * 3);

	for (size_t i = 0; i < vertexCount; ++i)
		vertexStartsFiber[i] = stream->readBool();
	vertexStartsFiber[vertexCount] = true;

	// TODO:
}

void BSplineHairShape::serialize(Stream *stream, InstanceManager *manager) const {
	Shape::serialize(stream, manager);
	Log(EError, "BSplineHairShape hasn't implemented serialize().");
	// TODO:
	// const std::vector<Point> &vertices = m_kdtree->getVertices();
	// const std::vector<bool> &vertexStartsFiber = m_kdtree->getStartFiber();

	// stream->writeFloat(m_kdtree->getRadius());
	// stream->writeSize(vertices.size());
	// stream->writeFloatArray((Float *)&vertices[0], vertices.size() * 3);
	// for (size_t i = 0; i < vertices.size(); ++i)
	//	stream->writeBool(vertexStartsFiber[i]);
}

bool BSplineHairShape::rayIntersect(const Ray &ray, Float mint,
									Float maxt, Float &t, void *temp) const {
	// return m_kdtree->rayIntersect(ray, mint, maxt, t, temp);
	static_assert(sizeof(RTCRayHit) <= MTS_KD_INTERSECTION_TEMP,
				  "RTCRayHit is too large for the temporary intersection buffer!");

	if (!m_embree_scene) {
		return false;
	}

	// Since we're not allocating, we can just use a local stack variable for the test.
	RTCRayHit rayhit;

	rayhit.ray.org_x = ray.o.x;
	rayhit.ray.org_y = ray.o.y;
	rayhit.ray.org_z = ray.o.z;
	rayhit.ray.dir_x = ray.d.x;
	rayhit.ray.dir_y = ray.d.y;
	rayhit.ray.dir_z = ray.d.z;
	rayhit.ray.tnear = mint;
	rayhit.ray.tfar = maxt;
	rayhit.ray.mask = -1;
	rayhit.ray.flags = 0;
	rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
	rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

	RTCIntersectArguments args;
	rtcInitIntersectArguments(&args);
	rtcIntersect1(m_embree_scene, &rayhit, &args);

	if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
		// A hit was found.
		t = rayhit.ray.tfar;
		if (temp != nullptr) {
			memcpy(temp, &rayhit, sizeof(RTCRayHit));
		}
		// Alternatively, a C++ style cast and assignment works too:
		// *static_cast<RTCRayHit*>(temp) = rayhit;
		return true;
	}

	return false;
}

bool BSplineHairShape::rayIntersect(const Ray &ray, Float mint, Float maxt) const {
	if (!m_embree_scene) {
		return false;
	}

	RTCRay embree_ray;
	embree_ray.org_x = ray.o.x;
	embree_ray.org_y = ray.o.y;
	embree_ray.org_z = ray.o.z;
	embree_ray.dir_x = ray.d.x;
	embree_ray.dir_y = ray.d.y;
	embree_ray.dir_z = ray.d.z;
	embree_ray.tnear = mint;
	embree_ray.tfar = maxt;
	embree_ray.mask = -1;
	embree_ray.flags = 0;

	RTCOccludedArguments args;
	rtcInitOccludedArguments(&args);
	rtcOccluded1(m_embree_scene, &embree_ray, &args);

	return embree_ray.tfar == -std::numeric_limits<float>::infinity();
}

void BSplineHairShape::fillIntersectionRecord(const Ray &ray,
											  const void *temp, Intersection &its) const {
	const RTCRayHit *rayhit = static_cast<const RTCRayHit *>(temp);

	its.shape = this;
	its.instance = this;
	its.time = ray.time;
	its.uv = Point2f(rayhit->hit.u, rayhit->hit.v);
	its.hasUVPartials = true;

	// Use rtcInterpolate to get the position and derivatives
	Point4f p;
	Vector4f dPdu, dPdv;
	rtcInterpolate1(m_embree_geom,
					rayhit->hit.primID,
					rayhit->hit.u,
					rayhit->hit.v,
					RTC_BUFFER_TYPE_VERTEX, 0, &p.x, &dPdu.x, &dPdv.x, 4);

	// Construct the geometric frame using data directly from the hit record
	// tangent:
	its.geoFrame.s = normalize(Vector3(dPdu.x, dPdu.y, dPdu.z));
	// normal:
	its.geoFrame.n = normalize(Normal(rayhit->hit.Ng_x, rayhit->hit.Ng_y, rayhit->hit.Ng_z));
	if (dot(its.geoFrame.n, ray.d) > 0) {
		its.geoFrame.n *= -1.0f;
	}
	// bitangent
	its.geoFrame.t = cross(its.geoFrame.n, its.geoFrame.s);

	its.shFrame = its.geoFrame;
	its.dpdu = its.geoFrame.s;
	its.dpdv = its.geoFrame.t;
	its.p = Point(p.x, p.y, p.z) + its.geoFrame.n * (p.w + Epsilon);

	its.splineT = rayhit->hit.u;
	its.primIndex = rayhit->hit.primID;
}

const std::vector<Point> &BSplineHairShape::getVertices() const {
	return m_original_vertices;
}

const std::vector<bool> &BSplineHairShape::getStartFiber() const {
	return m_original_vertex_starts_fiber;
}

const std::vector<uint32_t> &BSplineHairShape::getPrimitiveIndices() const {
	return m_embree_indices;
}

AABB BSplineHairShape::getAABB() const {
	if (!m_embree_scene) {
		return AABB(); // Return an empty bounding box if no geometry
	}

	RTCBounds bounds;
	rtcGetSceneBounds(m_embree_scene, &bounds);

	if (bounds.lower_x > bounds.upper_x) { // Check if bounds are invalid/empty
		return AABB();
	}

	return AABB(Point(bounds.lower_x, bounds.lower_y, bounds.lower_z),
				Point(bounds.upper_x, bounds.upper_y, bounds.upper_z));
}

// A "primitive" is one B-spline segment.
size_t BSplineHairShape::getPrimitiveCount() const {
	return m_embree_indices.size();
}

// For this shape, effective primitive count is the same.
size_t BSplineHairShape::getEffectivePrimitiveCount() const {
	return m_embree_indices.size();
}

bool BSplineHairShape::evaluatePrimitive(uint32_t primID, float u, CurvePointInfo &info) const {
	// 1. Input validation
	if (!m_embree_geom || primID >= m_embree_indices.size()) {
		// Geometry doesn't exist or primID is out of bounds
		return false;
	}

	// Clamp u to the valid [0, 1] range to be safe
	u = std::max(0.0f, std::min(1.0f, u));

	// 2. Prepare variables to receive interpolation results
	// Use Point4f for all to match the input buffer stride of 4 floats (x, y, z, radius)
	Point4f p_and_radius; // Embree will write x, y, z, radius here
	Point4f dPdU;		  // First derivative (tangent direction + dr/du)
	Point4f dPdV;		  // Second derivative/normal ( + dr/dv)

	// 3. Call rtcInterpolate1, which does all the hard work
	// We pass v=0 as requested. For a curve, v is not typically used.
	rtcInterpolate1(
		m_embree_geom,
		primID,
		u,
		0.0f, // v coordinate, set to 0
		RTC_BUFFER_TYPE_VERTEX, 0,
		&p_and_radius.x,
		&dPdU.x,
		&dPdV.x,
		4 // The stride of our vertex buffer is 4 floats (Point4f)
	);

	// 4. Populate the output structure
	info.p = Point(p_and_radius.x, p_and_radius.y, p_and_radius.z);
	info.radius = p_and_radius.w;

	// The first derivative gives the tangent direction (ignore w component).
	info.tangent = normalize(Vector(dPdU.x, dPdU.y, dPdU.z));

	// For curves, dPdV provides a vector perpendicular to the tangent,
	// which we can use as the bitangent (ignore w component).
	// info.bitangent = normalize(Vector(dPdV.x, dPdV.y, dPdV.z));

	// The normal can be found by taking the cross product.
	// This completes the right-handed Frenet frame (T, N, B).
	// info.normal = normalize(cross(info.bitangent, info.tangent));

	return true;
}

Point BSplineHairShape::getNextRayOrigin(const Intersection &its, const Vector &worldWo) const {
	if (dot(worldWo, its.geoFrame.n) < 0.0f) {
		return its.p - its.geoFrame.n.normalized() * (m_thickness + Epsilon);
	} else {
		return its.p;
	}
}

std::string BSplineHairShape::toString() const {
	std::ostringstream oss;
	// oss << "Hair[" << endl
	//	<< "   numVertices = " << m_kdtree->getVertexCount() << ","
	//	<< "   numSegments = " << m_kdtree->getSegmentCount() << ","
	//	<< "   numHairs = " << m_kdtree->getHairCount() << ","
	//	<< "   m_radius = " << m_kdtree->getRadius()
	//	<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S_API(BSplineHairShape, false, Shape, BSPLINEHAIR_API)
MTS_EXPORT_PLUGIN(BSplineHairShape, "Hair intersection shape");
MTS_NAMESPACE_END
