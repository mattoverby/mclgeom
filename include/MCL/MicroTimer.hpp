// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_MICROTIMER_HPP
#define MCL_MICROTIMER_HPP 1

#include <chrono>

namespace mcl
{

class MicroTimer
{
protected:
	//typedef std::chrono::high_resolution_clock C;
	typedef std::chrono::steady_clock C;
	typedef double T;

public:

	MicroTimer() : start_time(C::now()){}

	void reset() { start_time = C::now(); }

	// Return time elapsed in seconds
	T elapsed_s() const
	{
		curr_time = C::now();
		std::chrono::duration<T> durr = curr_time-start_time;
		return durr.count();
	}

	// Return time elapsed in milliseconds
	T elapsed_ms() const
	{
		curr_time = C::now();
		std::chrono::duration<T, std::milli> durr = curr_time-start_time;
		return durr.count();
	}

	// Return time elapsed in microseconds
	T elapsed_us() const
	{
		curr_time = C::now();
		std::chrono::duration<T, std::micro> durr = curr_time-start_time;
		return durr.count();
	}

private:
	std::chrono::time_point<C> start_time;
	mutable std::chrono::time_point<C> curr_time;

}; // end class MicroTimer

} // end namespace mcl

#endif
