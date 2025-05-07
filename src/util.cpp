#include <kklib/util.hpp>

void Timer::restart() { m_start = std::chrono::system_clock::now(); }

double Timer::duration()
{
    std::chrono::duration<double> diff = std::chrono::system_clock::now() - m_start;
    return diff.count();
}