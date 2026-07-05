#include <mitsuba/render/precomputer.h>

MTS_NAMESPACE_BEGIN

Precomputer::Precomputer(const Properties &props)
	: ConfigurableObject(props) {
	m_name = props.getID();
}

Precomputer::Precomputer(Stream *stream, InstanceManager *manager)
	: ConfigurableObject(stream, manager) {
}

MTS_IMPLEMENT_CLASS(Precomputer, true, ConfigurableObject)
MTS_NAMESPACE_END
