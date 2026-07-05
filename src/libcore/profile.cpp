#include <mitsuba/core/profiler.h>

MTS_NAMESPACE_BEGIN

MTS_EXPORT_CORE std::atomic<long long> IntegratorProfiler::m_accumulators[(size_t)ProfileCategory::COUNT];
MTS_EXPORT_CORE std::atomic<long long> IntegratorProfiler::m_executionCounters[(size_t)ProfileCategory::COUNT];
MTS_NAMESPACE_END