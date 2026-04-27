// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_MICROTIMER_HPP
#define MCL_GEOM_MICROTIMER_HPP 1

#include <chrono>

namespace mcl {

template<typename ClockType = std::chrono::steady_clock>
class MicroTimer
{
  protected:
    typedef double T;

  public:
    MicroTimer()
        : start_time(ClockType::now())
    {
    }

    void reset() { start_time = ClockType::now(); }

    // Return time elapsed in seconds
    T elapsed_s()
    {
        curr_time = ClockType::now();
        std::chrono::duration<T> durr = curr_time - start_time;
        return durr.count();
    }

    // Return time elapsed in milliseconds
    T elapsed_ms()
    {
        curr_time = ClockType::now();
        std::chrono::duration<T, std::milli> durr = curr_time - start_time;
        return durr.count();
    }

    // Return time elapsed in microseconds
    T elapsed_us()
    {
        curr_time = ClockType::now();
        std::chrono::duration<T, std::micro> durr = curr_time - start_time;
        return durr.count();
    }

  private:
    std::chrono::time_point<ClockType> start_time;
    std::chrono::time_point<ClockType> curr_time;

}; // end class MicroTimer

} // end namespace mcl

#endif
