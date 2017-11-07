#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <CL/cl.h>

#include "Dispatcher.hpp"
#include "ArgParser.hpp"
#include "constants.hpp"
#include "Mode.hpp"
#include "help.hpp"

std::string readFile(const char * const szFilename)
{
	std::ifstream in(szFilename, std::ios::in | std::ios::binary);
	std::ostringstream contents;
	contents << in.rdbuf();
	return contents.str();
}

std::vector<cl_device_id> getAllDevices(cl_device_type deviceType = CL_DEVICE_TYPE_GPU)
{
	std::vector<cl_device_id> vDevices;
	
	cl_uint platformIdCount = 0;
	clGetPlatformIDs (0, NULL, &platformIdCount);

	std::vector<cl_platform_id> platformIds (platformIdCount);
	clGetPlatformIDs (platformIdCount, platformIds.data (), NULL);
	
	for( auto it = platformIds.cbegin(); it != platformIds.cend(); ++it ) {
		cl_uint countDevice;
		clGetDeviceIDs(*it, deviceType, 0, NULL, &countDevice);
		
		std::vector<cl_device_id> deviceIds(countDevice);
		clGetDeviceIDs(*it, deviceType, countDevice, deviceIds.data(), &countDevice);
		
		std::copy( deviceIds.begin(), deviceIds.end(), std::back_inserter(vDevices) );
	}
	
	return vDevices;
}

template <typename T>
T clGetDeviceInfo(const cl_device_id & deviceId, const cl_device_info info) {
	T t;
	clGetDeviceInfo(deviceId, info, sizeof(t), &t, NULL);
	return t;
}

template <>
std::string clGetDeviceInfo(const cl_device_id & deviceId, const cl_device_info info) {
	size_t len;
	clGetDeviceInfo(deviceId, info, 0, NULL, &len);
	char * const szString = new char[len];
	clGetDeviceInfo(deviceId, info, len, szString, NULL);
	std::string r(szString);
	delete[] szString;
	return r;
}

template <typename T> bool printResult(const T & t, const cl_int & err) {
	std::cout << ((t == NULL) ? toString(err) : "OK") << std::endl;
	return t == NULL;
}

bool printResult(const cl_int err) {
	std::cout << ((err != CL_SUCCESS) ? toString(err) : "OK") << std::endl;
	return err != CL_SUCCESS;
}

int main(int argc, char * * argv) {
	try {
		ArgParser argp(argc, argv);
		bool bHelp = false;
		bool bModeBenchmark = false;
		bool bModeZeros = false;
		bool bModeLetters = false;
		bool bModeNumbers = false;
		std::string strModeLeading;
		std::string strModeMatching;
		bool bModeLeadingRange = false;
		bool bModeRange = false;
		int rangeMin = 0;
		int rangeMax = 0;
		std::vector<size_t> vDeviceSkipIndex;
		size_t worksizeLocal = 64;
		size_t worksizeMax = 32768;

		argp.addSwitch('h', "help", bHelp);
		argp.addSwitch('0', "benchmark", bModeBenchmark);
		argp.addSwitch('1', "zeros", bModeZeros);
		argp.addSwitch('2', "letters", bModeLetters);
		argp.addSwitch('3', "numbers", bModeNumbers);
		argp.addSwitch('4', "leading", strModeLeading);
		argp.addSwitch('5', "matching", strModeMatching);
		argp.addSwitch('6', "leading-range", bModeLeadingRange);
		argp.addSwitch('7', "range", bModeRange);
		argp.addSwitch('m', "min", rangeMin);
		argp.addSwitch('M', "max", rangeMax);
		argp.addMultiSwitch('s', "skip", vDeviceSkipIndex);
		argp.addSwitch('w', "work", worksizeLocal);
		argp.addSwitch('W', "workmax", worksizeMax);
		if (!argp.parse()) {
			std::cout << "error: bad arguments, try again :<" << std::endl;
			return 1;
		}

		if (bHelp) {
			std::cout << g_strHelp << std::endl;
			return 0;
		}

		Mode mode = Mode::benchmark();
		if (bModeBenchmark) {
			mode = Mode::benchmark();
		} else if (bModeZeros) {
			mode = Mode::zeros();
		} else if (bModeLetters) {
			mode = Mode::letters();
		} else if (bModeNumbers) {
			mode = Mode::numbers();
		} else if (!strModeLeading.empty()) {
			mode = Mode::leading(strModeLeading.front());
		} else if (!strModeMatching.empty()) {
			mode = Mode::matching(strModeMatching);
		} else if (bModeLeadingRange) {
			mode = Mode::leadingRange(rangeMin, rangeMax);
		} else if (bModeRange) {
			mode = Mode::range(rangeMin, rangeMax);
		} else {
			std::cout << g_strHelp << std::endl;
			return 0;
		}

		std::cout << "Mode: " << mode.name << std::endl;

		std::vector<cl_device_id> vDevicesFound = getAllDevices();
		cl_int errorCode;

		std::cout << "Devices:" << std::endl;
		for (int i = 0; i < vDevicesFound.size(); ++i) {
			cl_device_id & deviceId = vDevicesFound[i];

			const auto strName = clGetDeviceInfo<std::string>(deviceId, CL_DEVICE_NAME);
			const auto computeUnits = clGetDeviceInfo<cl_uint>(deviceId, CL_DEVICE_MAX_COMPUTE_UNITS);
			const auto globalMemSize = clGetDeviceInfo<cl_ulong>(deviceId, CL_DEVICE_GLOBAL_MEM_SIZE);

			std::cout << "\tGPU" << i << ": " << strName << ", " << globalMemSize << " bytes available, " << computeUnits << " compute units" << std::endl;
			std::cout << "\t\t" << "CL_DEVICE_MAX_CONSTANT_ARGS = " << clGetDeviceInfo<cl_uint>(deviceId, CL_DEVICE_MAX_CONSTANT_ARGS) << std::endl;
			std::cout << "\t\t" << "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE = " << clGetDeviceInfo<cl_ulong>(deviceId, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE) << std::endl;
			std::cout << "\t\t" << "CL_DEVICE_MAX_PARAMETER_SIZE = " << clGetDeviceInfo<size_t>(deviceId, CL_DEVICE_MAX_PARAMETER_SIZE) << std::endl;
			std::cout << std::endl;
		}

		// Device skip
		std::vector<cl_device_id> vDevices;
		for (size_t i = 0; i < vDevicesFound.size(); ++i) {
			if (std::find(vDeviceSkipIndex.begin(), vDeviceSkipIndex.end(), i) == vDeviceSkipIndex.end()) {
				vDevices.push_back(vDevicesFound[i]);
			}
		}

		if (vDevices.empty()) {
			return 1;
		}

		std::cout << "Initializing OpenCL..." << std::endl;
		std::cout << "\tCreating context..." << std::flush;
		auto clContext = clCreateContext( NULL, vDevices.size(), vDevices.data(), NULL, NULL, &errorCode);
		if (printResult(clContext, errorCode)) {
			return 1;
		}

		// Create a program from the kernel source
		const std::string strBignum = readFile("bignum.cl");
		const std::string strKeccak = readFile("keccak.cl");
		const std::string strVanity = readFile("profanity.cl");
		const char * szKernels[] = { strBignum.c_str(), strKeccak.c_str(), strVanity.c_str() };
	
		std::cout << "\tCompiling kernel..." << std::flush;
		cl_program clProgram = clCreateProgramWithSource(clContext, 3, szKernels, NULL, &errorCode);
		if (printResult(clProgram, errorCode)) {
			return 1;
		}
 
		// Build the program
		std::cout << "\tBuilding program..." << std::flush;
		const std::string strBuildOptions = "-D PROFANITY_PASSES=" + toString(PROFANITY_PASSES);
		if( printResult( clBuildProgram(clProgram, vDevices.size(), vDevices.data(), strBuildOptions.c_str(), NULL, NULL) ) ) {
#ifdef PROFANITY_DEBUG
			std::cout << std::endl;
			std::cout << "build log:" << std::endl;
		
			size_t sizeLog;
			clGetProgramBuildInfo(clProgram, vDevices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &sizeLog);
			char * const szLog = new char[sizeLog];
			clGetProgramBuildInfo(clProgram, vDevices[0], CL_PROGRAM_BUILD_LOG, sizeLog, szLog, NULL);
		
			std::cout << szLog << std::endl;
			delete[] szLog;
#endif
			return 1;
		}
		std::cout << std::endl;
 
		Dispatcher d(clContext, clProgram, mode, worksizeMax, 0);
		for (auto & i : vDevices) {
			d.addDevice(i, worksizeLocal);
		}
		
		d.run();
		clReleaseContext(clContext);
		return 0;
	} catch (std::runtime_error & e) {
		std::cout << "std::runtime_error - " << e.what() << std::endl;
	} catch (...) {
		std::cout << "unknown exception occured" << std::endl;
	}

	return 1;
}

