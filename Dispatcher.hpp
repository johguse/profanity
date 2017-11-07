#ifndef HPP_DISPATCHER
#define HPP_DISPATCHER

#include <stdexcept>
#include <fstream>
#include <CL/cl.h>
#include <string>
#include <chrono>
#include <mutex>
#include <list>
#include "types.hpp"
#include "CLMemory.hpp"
#include "Mode.hpp"

class Dispatcher {
	private:
		class OpenCLException : public std::runtime_error {
			public:
				OpenCLException(const std::string s, const cl_int res);

				static void throwIfError(const std::string s, const cl_int res);
		};

		struct SpeedSample {
			unsigned long long scanned;
			std::chrono::time_point<std::chrono::steady_clock> time;
		};
		
		struct Device {
			static cl_command_queue createQueue(cl_context & clContext, cl_device_id & clDeviceId);
			static cl_kernel createKernel(cl_program & clProgram, const std::string s);
			
			Device(Dispatcher & parent, cl_context & clContext, cl_program & clProgram, cl_device_id clDeviceId, const size_t worksizeLocal);
			~Device();
			
			Dispatcher & m_parent;
			
			cl_device_id m_clDeviceId;
			const size_t m_worksizeLocal;
			cl_uchar m_clScoreMax;
			cl_command_queue m_clQueue;
			
			cl_kernel m_kernelBegin;
			cl_kernel m_kernelInversePre;
			cl_kernel m_kernelInverse;
			cl_kernel m_kernelInversePost;
			cl_kernel m_kernelEnd;
			
			CLMemory<point> m_memPrecomp;
			CLMemory<point> m_memPoints1;
			CLMemory<point> m_memPoints2;
			CLMemory<bignum> m_memInverse;
			CLMemory<cl_uchar> m_memPass;
			
			CLMemory<result> m_memResult;
			
			// Data parameters used in some modes
			CLMemory<cl_uchar> m_memData1;
			CLMemory<cl_uchar> m_memData2;

			// Current seed for this device
			cl_ulong4 m_clSeed;
			
			std::list<SpeedSample> m_lSpeed;
			unsigned long long m_scanned;
		};
				
	public:
		Dispatcher(cl_context & clContext, cl_program & clProgram, const Mode mode, const size_t worksizeMax, const cl_uchar clScoreQuit = 0);
		~Dispatcher();
		
		void addDevice(cl_device_id clDeviceId, const size_t worksizeLocal);
		void run();
	
	private:
		void init(Device & d);
		void dispatch(Device & d);
		void enqueueKernel(cl_command_queue & clQueue, cl_kernel & clKernel, size_t worksizeGlobal, const size_t worksizeLocal);
		void handleResult(Device & d);
		void randomizeSeed(Device & d);
		
		void sampleAdd( std::list<SpeedSample> & lSpeed, const unsigned long long & scanned );
		unsigned int sampleSpeed( std::list<SpeedSample> & lSpeed );
		
		void onEvent(cl_event event, cl_int status, Device & d);
		
	private:
		static void CL_CALLBACK staticCallback(cl_event event, cl_int event_command_exec_status, void * user_data); 
		
		static std::string formatSpeed(float s);
		
	private: /* Instance variables */
		cl_context & m_clContext;
		cl_program & m_clProgram;
		const Mode m_mode;
		const size_t m_worksizeMax;
		cl_uchar m_clScoreMax;
		cl_uchar m_clScoreQuit;
		
		std::list<Device *> m_lDevices;
		
		cl_event m_eventFinished;

		// Run information
		std::mutex m_mutex;
		std::chrono::time_point<std::chrono::steady_clock> timeStart;
		std::list<SpeedSample> m_lSpeed;
		unsigned long long m_scanned;
		unsigned int m_countRunning;
		bool m_quit;
};

#endif /* HPP_DISPATCHER */