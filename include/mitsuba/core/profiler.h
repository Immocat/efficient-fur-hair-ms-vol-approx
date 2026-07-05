#pragma once

#include <atomic>
#include <boost/filesystem.hpp>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mitsuba/mitsuba.h>
#include <string>
#include <vector>

MTS_NAMESPACE_BEGIN

#define PROFILE_CATEGORIES(ACTION)                                                                                 \
	ACTION(Shading, "Full Shading")                                                                                \
	ACTION(SurfaceDirectLighting, "Surface Direct Shading")                                                        \
	ACTION(SurfaceSampling, "Surface Sampling")                                                                    \
	ACTION(DeltaTracking, "Delta Tracking")                                                                        \
	ACTION(VolumetricSingleScatterTransmittance, "Vol. SS Transmittance (RayMarching)")                            \
	ACTION(VolumetricMultiScatterScaleLoop, "Vol. MS Scaling Loop")                                                \
	ACTION(FullVolumetricShading, "Full Volumetric Shading")                                                       \
	ACTION(RenderBlockAccum, "Total RenderBlock() CPU Time")                                                       \
	ACTION(EvalEnvironment, "Time to eval environment in shading (only for those rays that did not hit geometry)") \
	ACTION(FirstIntersection, "Time to get the first intersection")                                                \
	ACTION(VolumetricLightSampling, "Time to sample a light source in the medium")

#define GENERATE_ENUM(name, string) name,
enum class ProfileCategory {
	PROFILE_CATEGORIES(GENERATE_ENUM)
		COUNT
};
#undef GENERATE_ENUM

#define GENERATE_STRING(name, string) string,
const std::vector<std::string> ProfileCategoryNames = {
	PROFILE_CATEGORIES(GENERATE_STRING)};
#undef GENERATE_STRING

class MTS_EXPORT_CORE IntegratorProfiler {
  public:
	static void reset() {
		for (size_t i = 0; i < (size_t)ProfileCategory::COUNT; ++i) {
			m_accumulators[i].store(0, std::memory_order_relaxed);
			m_executionCounters[i].store(0, std::memory_order_relaxed);
		}
	}

	static void add(ProfileCategory category, long long nanoseconds) {
		if (nanoseconds > 0) {
			m_accumulators[(size_t)category].fetch_add(nanoseconds, std::memory_order_relaxed);
			m_executionCounters[(size_t)category].fetch_add(1, std::memory_order_relaxed);
		}
	}

	static void report(float totalRenderTime, const fs::path &outputPath = fs::path()) {
		long long times_ns[(size_t)ProfileCategory::COUNT];
		long long counts[(size_t)ProfileCategory::COUNT];
		for (size_t i = 0; i < (size_t)ProfileCategory::COUNT; ++i) {
			times_ns[i] = m_accumulators[i].load();
			counts[i] = m_executionCounters[i].load(std::memory_order_relaxed);
		}

		using PC = ProfileCategory;

		long long shading_calls = counts[(size_t)PC::FirstIntersection];
		if (shading_calls == 0)
			shading_calls = 1;

		long long env_calls = counts[(size_t)PC::EvalEnvironment];
		long long hit_calls = shading_calls - env_calls;
		if (hit_calls < 0)
			hit_calls = 0;

		long long shading_other_ns = times_ns[(size_t)PC::Shading] - times_ns[(size_t)PC::FirstIntersection] - times_ns[(size_t)PC::EvalEnvironment] - times_ns[(size_t)PC::SurfaceDirectLighting] - times_ns[(size_t)PC::SurfaceSampling] - times_ns[(size_t)PC::FullVolumetricShading];
		if (shading_other_ns < 0)
			shading_other_ns = 0;

		long long vol_other_ns = times_ns[(size_t)PC::FullVolumetricShading] - times_ns[(size_t)PC::DeltaTracking]
								 - times_ns[(size_t)PC::VolumetricSingleScatterTransmittance] - times_ns[(size_t)PC::VolumetricMultiScatterScaleLoop];
		if (vol_other_ns < 0)
			vol_other_ns = 0;

		double total_cpu_s = (double)times_ns[(size_t)PC::RenderBlockAccum] / 1e9;
		double parallelism = (totalRenderTime > 0) ? (total_cpu_s / totalRenderTime) : 1.0;

		std::ofstream fileStream;
		bool useFile = !outputPath.empty();

		if (useFile) {
			try {
				fileStream.open(outputPath.string(), std::ios::out);
				if (!fileStream.is_open()) {
					SLog(EWarn, "Profiler: Could not open log file: %s", outputPath.string().c_str());
					useFile = false;
				}
			} catch (...) {
				SLog(EWarn, "Profiler: Exception while opening log file.");
				useFile = false;
			}
		}

		auto writeReport = [&](std::ostream &os) {
			os << "\n--- Integrator Performance Analysis ---\n";
			os << std::fixed << std::setprecision(2);
			os << "Timestamp              : " << std::time(nullptr) << " (Unix)\n";
			os << "Wall-clock Render Time : " << totalRenderTime << " s\n";
			os << "Total CPU Time         : " << total_cpu_s << " s\n";
			os << "Effective Parallelism  : " << parallelism << "x (Cores utilized)\n";
			os << "----------------------------------------------------------------------------------\n";
			os << std::left << std::setw(40) << "Category"
			   << std::right << std::setw(12) << "CPU Time"
			   << std::setw(12) << "Perc."
			   << std::setw(15) << "Est. Wall Time" << "\n";
			os << "----------------------------------------------------------------------------------\n";

			auto print_line = [&](const std::string &name, long long time_ns, long long parent_time_ns, const std::string &indent = "") {
				if (parent_time_ns == 0)
					return;

				double time_cpu_s = (double)time_ns / 1e9;
				double percentage = (double)time_ns / parent_time_ns * 100.0;
				double time_wall_s = time_cpu_s / parallelism;

				os << indent << std::left << std::setw(40 - indent.length()) << name
				   << ": " << std::right << std::setw(9) << time_cpu_s << " s"
				   << "  (" << std::setw(5) << percentage << "%)"
				   << "   ~" << std::setw(8) << time_wall_s << " s"
				   << std::endl;
			};

			print_line("Total Shading (Li)", times_ns[(size_t)PC::Shading], times_ns[(size_t)PC::RenderBlockAccum]);

			long long total_shading_time = times_ns[(size_t)PC::Shading];
			print_line("-> First Intersection", times_ns[(size_t)PC::FirstIntersection], total_shading_time, "  ");
			print_line("-> Environment Shading", times_ns[(size_t)PC::EvalEnvironment], total_shading_time, "  ");

			long long surface_and_vol_time = times_ns[(size_t)PC::SurfaceDirectLighting] + times_ns[(size_t)PC::SurfaceSampling] + times_ns[(size_t)PC::FullVolumetricShading];

			if (surface_and_vol_time > 0) {
				print_line("-> Surface & Volume (Total)", surface_and_vol_time, total_shading_time, "  ");
				print_line("  -> Surface Direct Light", times_ns[(size_t)PC::SurfaceDirectLighting], total_shading_time, "    ");
				print_line("  -> Surface BSDF Sample", times_ns[(size_t)PC::SurfaceSampling], total_shading_time, "    ");

				print_line("  -> Volumetric Shading", times_ns[(size_t)PC::FullVolumetricShading], total_shading_time, "    ");
				long long total_vol_time = times_ns[(size_t)PC::FullVolumetricShading];
				if (total_vol_time > 0) {
					print_line("    -> Delta Tracking", times_ns[(size_t)PC::DeltaTracking], total_vol_time, "      ");
					print_line("    -> Transmittance (RayMar)", times_ns[(size_t)PC::VolumetricSingleScatterTransmittance], total_vol_time, "      ");
					print_line("    -> Multi-Scatter Scale", times_ns[(size_t)PC::VolumetricMultiScatterScaleLoop], total_vol_time, "      ");
					print_line("    -> Volumetric Other", vol_other_ns, total_vol_time, "      ");
				}
			}

			print_line("-> Shading Logic & Other", shading_other_ns, total_shading_time, "  ");

			long long renderblock_other_ns = times_ns[(size_t)PC::RenderBlockAccum] - times_ns[(size_t)PC::Shading];
			if (renderblock_other_ns < 0)
				renderblock_other_ns = 0;

			os << "----------------------------------------------------------------------------------\n";
			print_line("RenderBlock Other", renderblock_other_ns, times_ns[(size_t)PC::RenderBlockAccum]);
			os << "----------------------------------------------------------------------------------\n\n";
		};

		writeReport(std::cout);

		if (useFile && fileStream.is_open()) {
			writeReport(fileStream);
			fileStream.close();
		}
	}

  private:
	static std::atomic<long long> m_accumulators[(size_t)ProfileCategory::COUNT];
	static std::atomic<long long> m_executionCounters[(size_t)ProfileCategory::COUNT];
};

class ScopeTimer {
  public:
	ScopeTimer(ProfileCategory category)
		: m_category(category), m_start(std::chrono::high_resolution_clock::now()) {
	}

	~ScopeTimer() {
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - m_start);
		IntegratorProfiler::add(m_category, duration.count());
	}

  private:
	ProfileCategory m_category;
	std::chrono::high_resolution_clock::time_point m_start;
};

class ManualTimer {
  public:
	explicit ManualTimer(ProfileCategory category)
		: m_category(category), m_isRunning(false) {
	}

	~ManualTimer() {
		if (m_isRunning) {
			SLog(EWarn, "ManualTimer for category '%s' was destroyed without calling stop()! Measurement lost.",
				 ProfileCategoryNames[(size_t)m_category].c_str());
		}
	}

	void start() {
		if (!m_isRunning) {
			m_start = std::chrono::high_resolution_clock::now();
			m_isRunning = true;
		}
	}

	void stop() {
		if (m_isRunning) {
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - m_start);
			IntegratorProfiler::add(m_category, duration.count());
			m_isRunning = false;
		}
	}

  private:
	ProfileCategory m_category;
	std::chrono::high_resolution_clock::time_point m_start;
	bool m_isRunning;
};

MTS_NAMESPACE_END
