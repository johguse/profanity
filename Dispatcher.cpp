#include "Dispatcher.hpp"

// Includes
#include <stdexcept>
#include <iostream>
#include <thread>
#include <sstream>
#include <iomanip>
#include <random>
#include <thread>
#include <algorithm>

#include "precomp.hpp"
#include "constants.hpp"

static std::string toHex(const uint8_t * const s, const size_t len) {
	std::string b("0123456789abcdef");
	std::string r;

	for (size_t i = 0; i < len; ++i) {
		const unsigned char h = s[i] / 16;
		const unsigned char l = s[i] % 16;

		r = r + b.substr(h, 1) + b.substr(l, 1);
	}

	return r;
}

static void printResult(cl_ulong4 seed, result r, const std::chrono::time_point<std::chrono::steady_clock> & timeStart) {
	// Time delta
	const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - timeStart).count();

	// Convert global finder id to bytes into part of seed
	cl_ulong seedPart = 0;
	for (int i = 0; i < PROFANITY_PASSES; ++i) {
		seedPart <<= 8;
		seedPart |= (r.foundId % 255) + 1;
		r.foundId /= 255;
	}

	// Mix seed.w with part of seed given above
	// Result: PPPSSSSS
	const auto bitSteps = 8 * (8 - PROFANITY_PASSES);
	const cl_ulong bitMask = (static_cast<cl_ulong>(1) << bitSteps) - 1;
	seed.s[3] = (seed.s[3] & bitMask) | (seedPart << bitSteps);

	// Format private key
	std::ostringstream ss;
	ss << std::hex << std::setfill('0');
	ss << std::setw(16) << seed.s[3] << std::setw(16) << seed.s[2] << std::setw(16) << seed.s[1] << std::setw(16) << seed.s[0];
	const std::string strPrivate = ss.str();

	// Format public key
	const std::string strPublic = toHex(r.foundHash, 20);

	// Print
	std::cout << "Time: " << std::setw(5) << seconds << "s Score: " << std::setw(2) << (int) r.foundScore << " Private: 0x" << strPrivate << " Public: 0x" << strPublic << std::endl;
}

Dispatcher::OpenCLException::OpenCLException(const std::string s, const cl_int res) :
	std::runtime_error( s + " (res = " + toString(res) + ")"),
	m_res(res)
{

}

void Dispatcher::OpenCLException::OpenCLException::throwIfError(const std::string s, const cl_int res) {
	if (res != CL_SUCCESS) {
		throw OpenCLException(s, res);
	}
}

cl_command_queue Dispatcher::Device::createQueue(cl_context & clContext, cl_device_id & clDeviceId) {
	const cl_command_queue ret = clCreateCommandQueueWithProperties(clContext, clDeviceId, NULL, NULL);
	return ret == NULL ? throw std::runtime_error("failed to create command queue") : ret;
}

cl_kernel Dispatcher::Device::createKernel(cl_program & clProgram, const std::string s) {
	cl_kernel ret  = clCreateKernel(clProgram, s.c_str(), NULL);
	return ret == NULL ? throw std::runtime_error("failed to create kernel \"" + s + "\"") : ret;
}

Dispatcher::Device::Device(Dispatcher & parent, cl_context & clContext, cl_program & clProgram, cl_device_id clDeviceId, const size_t worksizeLocal, const size_t index) :
	m_parent(parent),
	m_index(index),
	m_clDeviceId(clDeviceId),
	m_worksizeLocal(worksizeLocal),
	m_clScoreMax(0),
	m_clQueue(createQueue(clContext, clDeviceId) ),
	m_kernelBegin( createKernel(clProgram, "profanity_begin") ),
	m_kernelInversePre(createKernel(clProgram, "profanity_inverse_pre")),
	m_kernelInverse(createKernel(clProgram, "profanity_inverse_multiple")),
	m_kernelInversePost(createKernel(clProgram, "profanity_inverse_post")),
	m_kernelPass(createKernel(clProgram, "profanity_pass")),
	m_kernelEnd(createKernel(clProgram, "profanity_end")),
	m_memPrecomp(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, sizeof(g_precomp), g_precomp),
	m_memPoints(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, PROFANITY_MEM_SIZE, true),
	m_memPass(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, 1, true),
	m_memResult(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY, 40),
	m_memPointOffset(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, 1, true),
	m_memPointNextOffset(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, 1, true),
	m_memData1(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, 20),
	m_memData2(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, 20),
	m_speed(PROFANITY_SPEEDSAMPLES)
{

}

Dispatcher::Device::~Device() {

}

Dispatcher::Dispatcher(cl_context & clContext, cl_program & clProgram, const Mode mode, const size_t worksizeMax, const cl_uchar clScoreQuit)
	: m_clContext(clContext), m_clProgram(clProgram), m_mode(mode), m_worksizeMax(worksizeMax), m_clScoreMax(mode.score), m_clScoreQuit(clScoreQuit), m_eventFinished(NULL), m_countPrint(0) {

}

Dispatcher::~Dispatcher() {

}

void Dispatcher::addDevice(cl_device_id clDeviceId, const size_t worksizeLocal, const size_t index) {
	Device * pDevice = new Device(*this, m_clContext, m_clProgram, clDeviceId, worksizeLocal, index);
	m_lDevices.push_back(pDevice);
	init(*pDevice);
}

void Dispatcher::run() {
	m_eventFinished = clCreateUserEvent(m_clContext, NULL);
	if (m_eventFinished == NULL) {
		throw std::runtime_error("failed to create custom termination event");
	}

	m_quit = false;
	m_countRunning = m_lDevices.size();
	timeStart = std::chrono::steady_clock::now();

	for (auto it = m_lDevices.begin(); it != m_lDevices.end(); ++it) {
		dispatch(*(*it));
	}

	clWaitForEvents(1, &m_eventFinished);
	clReleaseEvent(m_eventFinished);
	m_eventFinished = NULL;
}

void Dispatcher::init(Device & d) {
	// Set mode data
	for (auto i = 0; i < 20; ++i) {
		d.m_memData1[i] = m_mode.data1[i];
		d.m_memData2[i] = m_mode.data2[i];
	}

	// Write precompute table and mode data
	d.m_memPrecomp.write(true);
	d.m_memData1.write(true);
	d.m_memData2.write(true);

	// Kernel arguments - profanity_begin
	d.m_memPrecomp.setKernelArg(d.m_kernelBegin, 0);
	d.m_memPoints.setKernelArg(d.m_kernelBegin, 1);
	d.m_memPointOffset.setKernelArg(d.m_kernelBegin, 2);
	d.m_memPointNextOffset.setKernelArg(d.m_kernelBegin, 3);
	d.m_memPass.setKernelArg(d.m_kernelBegin, 4);
	d.m_memResult.setKernelArg(d.m_kernelBegin, 5);
	/* seed set in dispatch() */

	// Kernel arguments - profanity_inverse_pre
	d.m_memPrecomp.setKernelArg(d.m_kernelInversePre, 0);
	d.m_memPoints.setKernelArg(d.m_kernelInversePre, 1);
	d.m_memPointOffset.setKernelArg(d.m_kernelInversePre, 2);
	d.m_memPointNextOffset.setKernelArg(d.m_kernelInversePre, 3);
	d.m_memPass.setKernelArg(d.m_kernelInversePre, 4);

	// Kernel arguments - profanity_inverse
	d.m_memPoints.setKernelArg(d.m_kernelInverse, 0);
	d.m_memPointNextOffset.setKernelArg(d.m_kernelInverse, 1);

	// Kernel arguments - profanity_inverse_post
	d.m_memPrecomp.setKernelArg(d.m_kernelInversePost, 0);
	d.m_memPoints.setKernelArg(d.m_kernelInversePost, 1);
	d.m_memPointOffset.setKernelArg(d.m_kernelInversePost, 2);
	d.m_memPointNextOffset.setKernelArg(d.m_kernelInversePost, 3);
	d.m_memPass.setKernelArg(d.m_kernelInversePost, 4);

	// Kernel arguments - profanity_pass
	d.m_memPass.setKernelArg(d.m_kernelPass, 0);
	d.m_memPointOffset.setKernelArg(d.m_kernelPass, 1);
	d.m_memPointNextOffset.setKernelArg(d.m_kernelPass, 2);

	// Kernel arguments - profanity_end
	d.m_memPoints.setKernelArg(d.m_kernelEnd, 0);
	d.m_memPointOffset.setKernelArg(d.m_kernelEnd, 1);
	d.m_memResult.setKernelArg(d.m_kernelEnd, 2);
	d.m_memData1.setKernelArg(d.m_kernelEnd, 3);
	d.m_memData2.setKernelArg(d.m_kernelEnd, 4);

	CLMemory<cl_uchar>::setKernelArg(d.m_kernelEnd, 5, d.m_clScoreMax);
	CLMemory<cl_uchar>::setKernelArg(d.m_kernelEnd, 6, m_mode.mode);
}

void Dispatcher::enqueueKernel(cl_command_queue & clQueue, cl_kernel & clKernel, size_t worksizeGlobal, const size_t worksizeLocal) {
	const size_t worksizeMax = m_worksizeMax;
	size_t worksizeOffset = 0;
	while (worksizeGlobal) {
		const size_t worksizeRun = std::min(worksizeGlobal, worksizeMax);
		const size_t * const pWorksizeLocal = (worksizeLocal == 0 ? NULL : &worksizeLocal);
		const auto res = clEnqueueNDRangeKernel(clQueue, clKernel, 1, &worksizeOffset, &worksizeRun, pWorksizeLocal, 0, NULL, NULL);
		OpenCLException::throwIfError("kernel queueing failed", res);

		worksizeGlobal -= worksizeRun;
		worksizeOffset += worksizeRun;
	}
}

void Dispatcher::enqueueKernelDevice(Device & d, cl_kernel & clKernel, size_t worksizeGlobal) {
	try {
		enqueueKernel(d.m_clQueue, clKernel, worksizeGlobal, d.m_worksizeLocal);
	}
	catch ( OpenCLException & e ) {
		// If local work size is invalid, abandon it and let implementation decide
		if ((e.m_res == CL_INVALID_WORK_GROUP_SIZE || e.m_res == CL_INVALID_WORK_ITEM_SIZE) && d.m_worksizeLocal != 0) {
			std::cout << std::endl << "warning: local work size abandoned on GPU" << d.m_index << std::endl;
			d.m_worksizeLocal = 0;
			enqueueKernel(d.m_clQueue, clKernel, worksizeGlobal, d.m_worksizeLocal);
		}
		else {
			throw;
		}
	}
}

void Dispatcher::dispatch(Device & d) {
	// Write new seed
	randomizeSeed(d);
	CLMemory<cl_ulong4>::setKernelArg(d.m_kernelBegin, 6, d.m_clSeed);

	enqueueKernelDevice(d, d.m_kernelBegin, 1);

	for (auto i = 1; i < PROFANITY_PASSES + 1; ++i) {
		enqueueKernelDevice(d, d.m_kernelInversePre, g_worksizes[i]);
		enqueueKernelDevice(d, d.m_kernelInverse, g_worksizes[i] / 255);
		enqueueKernelDevice(d, d.m_kernelInversePost, g_worksizes[i]);
		enqueueKernelDevice(d, d.m_kernelPass, 1);
	}

	enqueueKernelDevice(d, d.m_kernelEnd, g_worksizes[PROFANITY_PASSES]);

	cl_event event;
	d.m_memResult.read(false, &event);

	const auto res = clSetEventCallback(event, CL_COMPLETE, staticCallback, &d);
	OpenCLException::throwIfError("failed to set custom callback", res);
}

void Dispatcher::handleResult(Device & d) {
	//result & r = *d.m_memResult;
	for (auto i = 40; i > 0; --i) {
		result & r = d.m_memResult[i - 1];

		if (r.found > 0 && r.foundScore >= d.m_clScoreMax) {
			d.m_clScoreMax = r.foundScore;
			CLMemory<cl_uchar>::setKernelArg(d.m_kernelEnd, 5, d.m_clScoreMax);

			std::lock_guard<std::mutex> lock(m_mutex);
			if (r.foundScore >= m_clScoreMax) {
				m_clScoreMax = r.foundScore;

				if (m_clScoreQuit && r.foundScore >= m_clScoreQuit) {
					m_quit = true;
				}

				printResult(d.m_clSeed, r, timeStart);
			}

			break;
		}
	}
}

void Dispatcher::randomizeSeed(Device & d) {
	// Randomize private keys
	std::random_device rd;
	std::mt19937_64 eng(rd());
	std::uniform_int_distribution<cl_ulong> distr;
	d.m_clSeed.s[0] = distr(eng);
	d.m_clSeed.s[1] = distr(eng);
	d.m_clSeed.s[2] = distr(eng);
	d.m_clSeed.s[3] = distr(eng);
}

void Dispatcher::onEvent(cl_event event, cl_int status, Device & d) {
	if (status != CL_COMPLETE) {
		std::cout << "Dispatcher::onEvent - Got bad status: " << status << std::endl;
	}
	else {
		handleResult(d);

		bool bDispatch = true;
		{
			std::lock_guard<std::mutex> lock(m_mutex);
			d.m_speed.sample(PROFANITY_SIZE);
			printSpeed();

			if( m_quit ) {
				bDispatch = false;
				if(--m_countRunning == 0) {
					clSetUserEventStatus(m_eventFinished, CL_COMPLETE);
				}
			}
		}

		if (bDispatch) {
			dispatch(d);
		}
	}
}

// This is run when m_mutex is held.
void Dispatcher::printSpeed() {
	++m_countPrint;
	if( m_countPrint > m_lDevices.size() ) {
		std::string strGPUs;
		double speedTotal = 0;
		unsigned int i = 0;
		for (auto & e : m_lDevices) {
			const auto curSpeed = e->m_speed.getSpeed();
			speedTotal += curSpeed;
			strGPUs += " GPU" + toString(e->m_index) + ": " + formatSpeed(curSpeed);
			++i;
		}

		std::cout << "Total: " << formatSpeed(speedTotal) << " -" << strGPUs << "\r" << std::flush;
		m_countPrint = 0;
	}
}

void CL_CALLBACK Dispatcher::staticCallback(cl_event event, cl_int event_command_exec_status, void * user_data) {
	Device * const pDevice = static_cast<Device *>(user_data);
	pDevice->m_parent.onEvent(event, event_command_exec_status, *pDevice);
	clReleaseEvent(event);
}

std::string Dispatcher::formatSpeed(double f) {
	const std::string S = " KMGT";

	unsigned int index = 0;
	while (f > 1000.0f && index < S.size()) {
		f /= 1000.0f;
		++index;
	}

	std::ostringstream ss;
	ss << std::fixed << std::setprecision(3) << (double)f << " " << S[index] << "H/s";
	return ss.str();
}
