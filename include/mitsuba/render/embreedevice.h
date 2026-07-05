#pragma once
#if !defined(__EMBREEDEVICE_H_)
#define __EMBREEDEVICE_H_

#include <embree4/rtcore.h>
#include <iostream>

class EmbreeDevice {
  public:
	// Deleted copy constructor and assignment operator to prevent copying
	EmbreeDevice(const EmbreeDevice &) = delete;
	void operator=(const EmbreeDevice &) = delete;

	// Static method to get the single instance
	static EmbreeDevice &getInstance();

	// Method to get the shared device
	RTCDevice getDevice() {
		return m_device;
	}

  private:
	// Private constructor so it can only be created by getInstance()
	EmbreeDevice();

	// Private destructor to release the device on program exit
	~EmbreeDevice();

	// Static error handler function
	static void errorHandler(void *userPtr, RTCError error, const char *str);

	RTCDevice m_device = nullptr;
};

#endif