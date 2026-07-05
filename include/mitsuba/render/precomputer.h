#pragma once
#if !defined(__MITSUBA_RENDER_PRECOMPUTER_H_)
#define __MITSUBA_RENDER_PRECOMPUTER_H_
#include <mitsuba/core/cobject.h>
#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN
class MTS_EXPORT_RENDER Precomputer : public ConfigurableObject {
  public:
	/**
	 * \brief Phase 1: Pre-Build
	 *
	 * This method is called BEFORE the scene's acceleration structure (KD-tree)
	 * is constructed.
	 *
	 * \note Use this to add geometry, modify meshes, or change material assignments.
	 * \warning You cannot cast rays or query the scene for intersections here.
	 */
	virtual void runPreBuild(Scene *scene) {
	}

	/**
	 * \brief Phase 2: Post-Build
	 *
	 * This method is called AFTER the scene's acceleration structure is fully
	 * initialized.
	 *
	 * \note Use this for tasks requiring ray tracing, such as irradiance caching,
	 * baking textures, or voxelization.
	 */
	virtual void runPostBuild(Scene *scene) {
	}

	Precomputer(const Properties &props);

	Precomputer(Stream *stream, InstanceManager *manager);

	MTS_DECLARE_CLASS()
  protected:
	virtual ~Precomputer() {
	}

	std::string m_name;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_PRECOMPUTER_H_ */