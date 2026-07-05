#include <mitsuba/render/embreedevice.h>
// Static method to get the single instance
EmbreeDevice &EmbreeDevice::getInstance() {
	static EmbreeDevice instance;
	return instance;
}

// Private constructor
EmbreeDevice::EmbreeDevice() {
	m_device = rtcNewDevice(NULL);
	if (!m_device) {
		// Handle error, maybe throw an exception
		std::cerr << "Failed to create Embree device!" << std::endl;
	}

	// Set a global error handler
	rtcSetDeviceErrorFunction(m_device, errorHandler, nullptr);
}

// Private destructor
EmbreeDevice::~EmbreeDevice() {
	if (m_device) {
		rtcReleaseDevice(m_device);
	}
}

// Static error handler function
void EmbreeDevice::errorHandler(void *userPtr, RTCError error, const char *str) {
	std::cerr << "Embree Error: " << str << std::endl;
}