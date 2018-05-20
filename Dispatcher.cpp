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

static void printResult(cl_ulong4 seed, cl_ulong round, result r, cl_uchar score, const std::chrono::time_point<std::chrono::steady_clock> & timeStart) {
	// Time delta
	const auto seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - timeStart).count();

	// Format private key
	cl_ulong carry = 0;
	cl_ulong4 seedRes;

	seedRes.s[0] = seed.s[0] + round; carry = seedRes.s[0] < round;
	seedRes.s[1] = seed.s[1] + carry; carry = !seedRes.s[1];
	seedRes.s[2] = seed.s[2] + carry; carry = !seedRes.s[2];
	seedRes.s[3] = seed.s[3] + carry + r.foundId;

	std::ostringstream ss;
	ss << std::hex << std::setfill('0');
	ss << std::setw(16) << seedRes.s[3] << std::setw(16) << seedRes.s[2] << std::setw(16) << seedRes.s[1] << std::setw(16) << seedRes.s[0];
	const std::string strPrivate = ss.str();

	// Format public key
	const std::string strPublic = toHex(r.foundHash, 20);

	// Print
	std::cout << "  Time: " << std::setw(5) << seconds << "s Score: " << std::setw(2) << (int) score << " Private: 0x" << strPrivate << " Public: 0x" << strPublic << std::endl;
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

cl_ulong4 Dispatcher::Device::createSeed() {
	// Randomize private keys
	std::random_device rd;
	std::mt19937_64 eng(rd());
	std::uniform_int_distribution<cl_ulong> distr;

	cl_ulong4 r;
	r.s[0] = distr(eng);
	r.s[1] = distr(eng);
	r.s[2] = distr(eng);
	r.s[3] = distr(eng);
	return r;
}

Dispatcher::Device::Device(Dispatcher & parent, cl_context & clContext, cl_program & clProgram, cl_device_id clDeviceId, const size_t worksizeLocal, const size_t size, const size_t index, const Mode & mode) :
	m_parent(parent),
	m_index(index),
	m_clDeviceId(clDeviceId),
	m_worksizeLocal(worksizeLocal),
	m_clScoreMax(0),
	m_clQueue(createQueue(clContext, clDeviceId) ),
	m_kernelBegin( createKernel(clProgram, "profanity_begin") ),
	m_kernelInverse(createKernel(clProgram, "profanity_inverse_multiple")),
	m_kernelInversePost(createKernel(clProgram, "profanity_inverse_post")),
	m_kernelEnd(createKernel(clProgram, "profanity_end")),
	m_kernelScore(createKernel(clProgram, mode.kernel)),
	m_memPrecomp(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, sizeof(g_precomp), g_precomp),
	m_memPoints(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, size, true),
	m_memInverse(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS, size, true),
	m_memResult(clContext, m_clQueue, CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY, 40),
	m_memData1(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, 20),
	m_memData2(clContext, m_clQueue, CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY, 20),
	m_clSeed(createSeed()),
	m_round(0),
	m_speed(PROFANITY_SPEEDSAMPLES),
	m_sizeInitialized(0),
	m_eventFinished(NULL)
{

}

Dispatcher::Device::~Device() {

}

Dispatcher::Dispatcher(cl_context & clContext, cl_program & clProgram, const Mode mode, const size_t worksizeMax, const size_t inverseSize, const size_t inverseMultiple, const cl_uchar clScoreQuit)
	: m_clContext(clContext), m_clProgram(clProgram), m_mode(mode), m_worksizeMax(worksizeMax), m_inverseSize(inverseSize), m_size(inverseSize*inverseMultiple), m_clScoreMax(mode.score), m_clScoreQuit(clScoreQuit), m_eventFinished(NULL), m_countPrint(0) {

}

Dispatcher::~Dispatcher() {

}

void Dispatcher::addDevice(cl_device_id clDeviceId, const size_t worksizeLocal, const size_t index) {
	Device * pDevice = new Device(*this, m_clContext, m_clProgram, clDeviceId, worksizeLocal, m_size, index, m_mode);
	m_vDevices.push_back(pDevice);
}

void Dispatcher::run() {
	m_eventFinished = clCreateUserEvent(m_clContext, NULL);
	init();

	m_quit = false;
	m_countRunning = m_vDevices.size();
	timeStart = std::chrono::steady_clock::now();

	std::cout << "Running..." << std::endl;
	std::cout << "  Always verify that a private key generated by this program corresponds to the" << std::endl;
	std::cout << "  public key printed by importing it to a wallet of your choice. This program" << std::endl;
	std::cout << "  like any software might contain bugs and it does by design cut corners to" << std::endl;
	std::cout << "  improve overall performance." << std::endl;
	std::cout << std::endl;

	for (auto it = m_vDevices.begin(); it != m_vDevices.end(); ++it) {
		dispatch(*(*it));
	}

	clWaitForEvents(1, &m_eventFinished);
	clReleaseEvent(m_eventFinished);
	m_eventFinished = NULL;
}

void Dispatcher::init() {
	std::cout << "Initializing devices..." << std::endl;
	std::cout << "  This can take a minute or two. The number of objects initialized on each" << std::endl;
	std::cout << "  device is equal to inverse-size * inverse-multiple. To lower" << std::endl;
	std::cout << "  initialization time (and memory footprint) I suggest lowering the" << std::endl;
	std::cout << "  inverse-multiple first. You can do this via the -I switch. Do note that" << std::endl;
	std::cout << "  this might negatively impact your performance." << std::endl;
	std::cout << std::endl;

	const auto deviceCount = m_vDevices.size();

	cl_event * const pInitEvents = new cl_event[deviceCount];

	for (size_t i = 0; i < deviceCount; ++i) {
		pInitEvents[i] = clCreateUserEvent(m_clContext, NULL);
		m_vDevices[i]->m_eventFinished = pInitEvents[i];
		initBegin(*m_vDevices[i]);
	}

	clWaitForEvents(deviceCount, pInitEvents);
	for (size_t i = 0; i < deviceCount; ++i) {
		m_vDevices[i]->m_eventFinished = NULL;
		clReleaseEvent(pInitEvents[i]);
	}

	delete[] pInitEvents;

	std::cout << std::endl;
}

void Dispatcher::initBegin(Device & d) {
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
	d.m_memResult.setKernelArg(d.m_kernelBegin, 2);
	CLMemory<cl_ulong4>::setKernelArg(d.m_kernelBegin, 3, d.m_clSeed);

	// Kernel arguments - profanity_inverse
	d.m_memPoints.setKernelArg(d.m_kernelInverse, 0);
	d.m_memInverse.setKernelArg(d.m_kernelInverse, 1);

	// Kernel arguments - profanity_inverse_post
	d.m_memPoints.setKernelArg(d.m_kernelInversePost, 0);
	d.m_memInverse.setKernelArg(d.m_kernelInversePost, 1);

	// Kernel arguments - profanity_end
	d.m_memPoints.setKernelArg(d.m_kernelEnd, 0);
	d.m_memInverse.setKernelArg(d.m_kernelEnd, 1);

	// Kernel arguments - profanity_score_???
	d.m_memInverse.setKernelArg(d.m_kernelScore, 0);
	d.m_memResult.setKernelArg(d.m_kernelScore, 1);
	d.m_memData1.setKernelArg(d.m_kernelScore, 2);
	d.m_memData2.setKernelArg(d.m_kernelScore, 3);

	CLMemory<cl_uchar>::setKernelArg(d.m_kernelScore, 4, d.m_clScoreMax); // Updated in handleResult()

	// Seed device
	initContinue(d);
}

void Dispatcher::initContinue(Device & d) {
	size_t sizeLeft = m_size - d.m_sizeInitialized;

	if (sizeLeft) {
		cl_event event;
		const size_t sizeRun = std::min(sizeLeft, m_worksizeMax);
		const auto resEnqueue = clEnqueueNDRangeKernel(d.m_clQueue, d.m_kernelBegin, 1, &d.m_sizeInitialized, &sizeRun, NULL, 0, NULL, &event);
		OpenCLException::throwIfError("kernel queueing failed during initilization", resEnqueue);

		const auto resCallback = clSetEventCallback(event, CL_COMPLETE, staticCallback, &d);
		OpenCLException::throwIfError("failed to set custom callback during initialization", resCallback);

		d.m_sizeInitialized += sizeRun;
	} else {
		// Printing one whole string at once helps in avoiding garbled output when executed in parallell
		const std::string strOutput = "  GPU" + toString(d.m_index) + " initialized";
		std::cout << strOutput << std::endl;
		clSetUserEventStatus(d.m_eventFinished, CL_COMPLETE);
	}
}

void Dispatcher::enqueueKernel(cl_command_queue & clQueue, cl_kernel & clKernel, size_t worksizeGlobal, const size_t worksizeLocal, const bool bOneAtATime = false) {
	const size_t worksizeMax = m_worksizeMax;
	size_t worksizeOffset = 0;
	cl_event clEvent;
	while (worksizeGlobal) {
		const size_t worksizeRun = std::min(worksizeGlobal, worksizeMax);
		const size_t * const pWorksizeLocal = (worksizeLocal == 0 ? NULL : &worksizeLocal);
		const auto res = clEnqueueNDRangeKernel(clQueue, clKernel, 1, &worksizeOffset, &worksizeRun, pWorksizeLocal, 0, NULL, bOneAtATime ? &clEvent : NULL);
		OpenCLException::throwIfError("kernel queueing failed", res);

		// Queueing lots of work exhausted resources on my GTX 1070 during initialization. I don't really know why. Correlated with worksizeMax.
		if (bOneAtATime) {
			clWaitForEvents(1, &clEvent);
			clReleaseEvent(clEvent);
			clEvent = NULL;
		}

		worksizeGlobal -= worksizeRun;
		worksizeOffset += worksizeRun;
	}
}

void Dispatcher::enqueueKernelDevice(Device & d, cl_kernel & clKernel, size_t worksizeGlobal, const bool bOneAtATime = false) {
	try {
		enqueueKernel(d.m_clQueue, clKernel, worksizeGlobal, d.m_worksizeLocal, bOneAtATime);
	}
	catch ( OpenCLException & e ) {
		// If local work size is invalid, abandon it and let implementation decide
		if ((e.m_res == CL_INVALID_WORK_GROUP_SIZE || e.m_res == CL_INVALID_WORK_ITEM_SIZE) && d.m_worksizeLocal != 0) {
			std::cout << std::endl << "warning: local work size abandoned on GPU" << d.m_index << std::endl;
			d.m_worksizeLocal = 0;
			enqueueKernel(d.m_clQueue, clKernel, worksizeGlobal, d.m_worksizeLocal, bOneAtATime);
		}
		else {
			throw;
		}
	}
}

void Dispatcher::dispatch(Device & d) {
	enqueueKernelDevice(d, d.m_kernelInverse, m_size / m_inverseSize);
	enqueueKernelDevice(d, d.m_kernelInversePost, m_size);
	enqueueKernelDevice(d, d.m_kernelEnd, m_size);
	enqueueKernelDevice(d, d.m_kernelScore, m_size);

	cl_event event;
	d.m_memResult.read(false, &event);

	const auto res = clSetEventCallback(event, CL_COMPLETE, staticCallback, &d);
	OpenCLException::throwIfError("failed to set custom callback", res);
}

void Dispatcher::handleResult(Device & d) {
	++d.m_round;

	for (auto i = 40 - 1; i > m_clScoreMax; --i) {
		result & r = d.m_memResult[i];

		if (r.found > 0 && i >= d.m_clScoreMax) {
			d.m_clScoreMax = i;
			CLMemory<cl_uchar>::setKernelArg(d.m_kernelScore, 4, d.m_clScoreMax);

			std::lock_guard<std::mutex> lock(m_mutex);
			if (i >= m_clScoreMax) {
				m_clScoreMax = i;

				if (m_clScoreQuit && i >= m_clScoreQuit) {
					m_quit = true;
				}

				printResult(d.m_clSeed, d.m_round, r, i, timeStart);
			}

			break;
		}
	}
}

void Dispatcher::onEvent(cl_event event, cl_int status, Device & d) {
	if (status != CL_COMPLETE) {
		std::cout << "Dispatcher::onEvent - Got bad status: " << status << std::endl;
	}
	else if (d.m_eventFinished != NULL) {
		initContinue(d);
	} else {
		handleResult(d);

		bool bDispatch = true;
		{
			std::lock_guard<std::mutex> lock(m_mutex);
			d.m_speed.sample(m_size);
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
	if( m_countPrint > m_vDevices.size() ) {
		std::string strGPUs;
		double speedTotal = 0;
		unsigned int i = 0;
		for (auto & e : m_vDevices) {
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
