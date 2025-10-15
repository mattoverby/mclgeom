// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SMOOTHMOVE_HPP
#define MCL_SMOOTHMOVE_HPP 1

#include <Eigen/Dense>

namespace mcl {

// Useful for moving control points smoothly
// total_elapsed_dt = system time elapsed sum
// start_dt = time to start the movement
// end_dt = time the movement should stop

template<typename T>
inline Eigen::Vector3<T>
smooth_move(T total_elapsed_dt, T start_dt, T end_dt, const Eigen::Vector3<T>& start, const Eigen::Vector3<T>& end)
{
    if (total_elapsed_dt < start_dt) {
        return start;
    }
    T t_ratio = (total_elapsed_dt - start_dt) / (end_dt - start_dt);
    if (t_ratio > 1.0) {
        return end;
    }
    Eigen::Vector3<T> displacement = end - start;
    return (start + (3.0 * t_ratio * t_ratio - 2.0 * t_ratio * t_ratio * t_ratio) * displacement);
}

template<typename T>
inline Eigen::Vector3<T>
linear_move(T total_elapsed_dt, T start_dt, T end_dt, const Eigen::Vector3<T>& start, const Eigen::Vector3<T>& end)
{
    if (total_elapsed_dt < start_dt) {
        return start;
    }
    T t_ratio = (total_elapsed_dt - start_dt) / (end_dt - start_dt);
    if (t_ratio > 1.0) {
        return end;
    }
    Eigen::Vector3<T> displacement = end - start;
    return (start + displacement);
}

} // namespace mcl

#endif