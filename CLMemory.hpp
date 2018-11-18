#ifndef HPP_CLMEMORY
#define HPP_CLMEMORY

#include "lexical_cast.hpp"

template<typename T> class CLMemory {
	public:
		CLMemory(cl_context & clContext, cl_command_queue & clQueue, const cl_mem_flags flags, const size_t size, T * const pData)
		 : m_clQueue(clQueue), m_bFree(false), m_size(size), m_pData(pData) {
			 m_clMem = clCreateBuffer(clContext, flags, m_size, NULL, NULL);
		}

		CLMemory(cl_context & clContext, cl_command_queue & clQueue, const cl_mem_flags flags, const size_t count, const bool noAllocation = false)
		 : m_clQueue(clQueue), m_bFree(true), m_size(sizeof(T) * count), m_pData(noAllocation ? NULL : new T[count]) {
			m_clMem = clCreateBuffer(clContext, flags, m_size, NULL, NULL);
		}

		~CLMemory() {
			if(m_bFree) {
				delete [] m_pData;
			}
		}

		static void setKernelArg(cl_kernel & clKernel, const cl_uint arg_index, const T & t) {
			const cl_int ret = clSetKernelArg(clKernel, arg_index, sizeof(T), (void *) &t);
			if (ret != CL_SUCCESS) {
				throw std::runtime_error("clSetKernelArg failed - " + toString(arg_index) + " - " + toString(ret));
			}
		}

		void setKernelArg(cl_kernel & clKernel, const cl_uint arg_index) const {
			const cl_int ret = clSetKernelArg(clKernel, arg_index, sizeof(cl_mem), (void *) &m_clMem );
			if( ret != CL_SUCCESS ) {
				throw std::runtime_error("clSetKernelArg failed - " + toString(arg_index) + " - " + toString(ret));
			}
		}

		void read(const bool bBlock, cl_event * pEvent = NULL) const {
			const cl_bool block = bBlock ? CL_TRUE : CL_FALSE;
			auto res = clEnqueueReadBuffer(m_clQueue, m_clMem, block, 0, m_size, m_pData, 0, NULL, pEvent);
			if(res != CL_SUCCESS) {
				throw std::runtime_error("clEnqueueReadBuffer failed - " + toString(res));
			}
		}

		void write(const bool bBlock) const {
			const cl_bool block = bBlock ? CL_TRUE : CL_FALSE;
			auto res = clEnqueueWriteBuffer(m_clQueue, m_clMem, block, 0, m_size, m_pData, 0, NULL, NULL);
			if( res != CL_SUCCESS ) {
				throw std::runtime_error("clEnqueueWriteBuffer failed - " + toString(res));
			}
		}

		T * const & data() const {
			return m_pData;
		}

		T & operator[] (int x) const {
			return m_pData[x];
		}

		T * operator->() const {
			return m_pData;
		}

		T & operator*() const {
			return *m_pData;
		}

		const size_t & size() const {
			return m_size;
		}

	private:
		const cl_command_queue m_clQueue;
		const bool m_bFree;
		const size_t m_size;

		T * const m_pData;
		cl_mem m_clMem;
};

#endif /* HPP_CLMEMORY */