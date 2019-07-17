#include "SpeedSample.hpp"

SpeedSample::SpeedSample(const size_t length) :
	m_length(length),
	m_lastTime(now())
{

}

SpeedSample::~SpeedSample() {

}

double SpeedSample::getSpeed() const {
	auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(now() - m_lastTime).count();
	if (delta > 5000) {
		return 0;
	} else {
		double speed = 0;
		for (auto & v : m_lSpeeds) {
			speed += v / m_lSpeeds.size();
		}

		return speed;
	}
}

void SpeedSample::sample(const double V) {
	const timepoint newTime = now();
	auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(newTime - m_lastTime).count();
	m_lSpeeds.push_back((1000 * V) / delta);
	m_lastTime = newTime;
	if (m_lSpeeds.size() > m_length) {
		m_lSpeeds.pop_front();
	}
}

SpeedSample::timepoint SpeedSample::now() {
	return std::chrono::steady_clock::now();
}
