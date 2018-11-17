#ifndef HPP_SPEEDSAMPLE
#define HPP_SPEEDSAMPLE
#include <chrono>
#include <list>

class SpeedSample {
	private:
		typedef std::chrono::time_point<std::chrono::steady_clock> timepoint;

	public:
		SpeedSample(const size_t length);
		~SpeedSample();

		double getSpeed() const;
		void sample(const double V);

	private:
		static timepoint now();

	private:
		const size_t m_length;
		timepoint m_lastTime;
		std::list<double> m_lSpeeds;
};

#endif /* HPP_SPEEDSAMPLE */