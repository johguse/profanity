#ifndef HPP_MODE
#define HPP_MODE

#include <string>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

class Mode {
	private:
		Mode();

	public:
		static Mode matching(const std::string strHex);
		static Mode range(const cl_uchar min, const cl_uchar max);
		static Mode leading(const char charLeading);
		static Mode leadingRange(const cl_uchar min, const cl_uchar max);
		static Mode mirror();

		static Mode benchmark();
		static Mode zeros();
		static Mode letters();
		static Mode numbers();

		std::string name;
		std::string kernel;

		cl_uchar data1[20];
		cl_uchar data2[20];
		cl_uchar score;
};

#endif /* HPP_MODE */