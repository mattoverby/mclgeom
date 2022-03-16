// Copyright Matt Overby 2022.
// Distributed under the MIT License.

#ifndef MCL_LOGGER_HPP
#define MCL_LOGGER_HPP 1

#include "MCL/MicroTimer.hpp"
#include <mutex>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iomanip>
#include <fstream>

namespace mcl
{

// This is a very simply logging tool for recording run time on a per-frame basis.
// Macros support this functionality, and the underlying class is a singleton.
//
// Examples:
//
// void Solver::init()
// {
//   mcl::ResetLog();
//   mcl::TimedScope("Solver::init");
//   init_the_solver();
// } // run time of the function is recorded
//
// void Solver::solve()
// {
//   mcl::IncrementFrame();
//   solve_the_frame();
//   mcl::StopFrame(); (optional, stops timer and records run time)
//   lots_of_expensive_logging();
// } // run time is recorded here if StopFrame() not called.
//
// If not all computation is captured in a single scope (like above):
//   mcl::AddRuntime(elapsed_ms); // adds run time to current frame
//   mcl::AddFrameRuntime(elapsed_ms, frame); // adds time to specific frame
//
// You can also time a specific scope/function:
// {
//   mcl::TimedScope("MyClass::MyFunction");
//   ...
// } // timer stops, and function run time is recorded
//
// If timers are inside of a loop, the total is summed:
// for (int i=0; i<10; ++i)
// {
//   mcl::TimedScope("MyClass::MyLoop");
//   ...
// } // total time for all 10 iterations is accumulated
//
// The timings can be written to CSV file:
//   mcl::WriteLogCSV("mysolver_log.csv");
//
class FunctionTimerWrapper;
#define ResetLog() FunctionTimerWrapper::reset_all()
#define IncrementFrame() FunctionTimerWrapper::increment_frame()
#define StopFrame() FunctionTimerWrapper::stop_frame();
#define	AddRuntime(ms) FunctionTimerWrapper::add_runtime(ms, -2)
#define	AddFrameRuntime(ms, frame) FunctionTimerWrapper::add_runtime(ms, frame)
#define TimedScope(f) FunctionTimerWrapper ftWrapper_(f)
#define WriteLogCSV(csv) FunctionTimerWrapper::write_csv(csv)

// Singleton class for Logging
// e.g., mcl::Logger &log = mcl::Logger::get();
//
class Logger
{
protected:
	Logger() {}
	virtual ~Logger() {}

public:
	Logger(Logger const &) = delete;
	Logger operator=(Logger const &) = delete;

	static Logger& get()
	{
		static Logger instance;
		return instance;
	}

	// Clears all logged data
	void clear()
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		curr_frame = -1;
		init_data = PerFrameData();
		frame_data.clear();
	}

	int frame_number() const { return curr_frame; }

	// Increases the frame count and begin recording at new frame.
	void increment_frame()
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		curr_frame++;
		while ((int)frame_data.size() < curr_frame+1) { frame_data.emplace_back(); }
		if (curr_frame >= 0) { frame_data[curr_frame].start_timer("elapsed_ms"); }
	}

	void stop_frame()
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		while ((int)frame_data.size() < curr_frame+1) { frame_data.emplace_back(); }
		if (curr_frame >= 0) { frame_data[curr_frame].stop_all_timers(); }
	}

	// Variable used for keeping track of the per-frame runtime. If -2, uses current frame.
	void add_runtime_s(double s, int frame = -2)
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		if (frame == -2) { frame = curr_frame; }
		if (curr_frame == -1) { init_data.add_time("elapsed_ms", s); }
		else if (curr_frame >= 0)
		{
			while ((int)frame_data.size() < frame+1) { frame_data.emplace_back(); }
			frame_data[curr_frame].add_time("elapsed_ms", s);
		}
	}

	void start_timer(const std::string &l)
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		if (curr_frame == -1) { init_data.start_timer(l); }
		else if (curr_frame >= 0)
		{
			while ((int)frame_data.size() < curr_frame+1) { frame_data.emplace_back(); }
			frame_data[curr_frame].start_timer(l);
		}
	}

	void stop_timer(const std::string &l)
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		if (curr_frame == -1) { init_data.stop_timer(l); }
		else if (curr_frame >= 0)
		{
			while ((int)frame_data.size() < curr_frame+1) { frame_data.emplace_back(); }
			frame_data[curr_frame].stop_timer(l);
		}
	}

	// Returns per-frame (excluding init) timings in millesconds for every timer.
	// If a timing wasn't recorded for a given frame, it is set to zero.
	// e.g., double time_ms = timings[label][frame_number].
	void get_timings(std::unordered_map<std::string, std::vector<double> > &timings_ms)
	{
		std::lock_guard<std::mutex> guard(write_mutex);
		int num_data = frame_data.size();

		// Get list of all labels. This ensures there is a value at every frame.
		// Also insert all labels into init_data (as zero if not found).
		typedef std::unordered_map<std::string, int>::iterator LabelIt;
		std::unordered_map<std::string, int> labels;
		for (int i=-1; i<num_data; ++i)
		{
			PerFrameData &data = (i==-1) ? init_data : frame_data[i];
			data.stop_all_timers();
			for (DoubleIt it = data.timings_ms.begin(); it != data.timings_ms.end(); ++it)
			{
				labels.emplace(it->first, 0);
				init_data.timings_ms.emplace(it->first, 0);
			}
		}

		// Loop through per-frame data and insert a zero if the timer is not present in the frame.
		for (int i=0; i<num_data; ++i)
		{
			for (LabelIt it = labels.begin(); it != labels.end(); ++it)
			{
				// Insert zeros at empty spaces up to current frame
				std::vector<double> &per_frame_ms = timings_ms.emplace(it->first, std::vector<double>()).first->second;
				while ((int)per_frame_ms.size() < i+1) { per_frame_ms.emplace_back(0); }
				// Set elapsed if present in the frame data
				DoubleIt it_ms = frame_data[i].timings_ms.find(it->first);
				if (it_ms != frame_data[i].timings_ms.end()) { per_frame_ms[i] = it_ms->second; }
			}
		}
	}

	// Saves logged data to a CSV.
	// Total elapsed time is summed per frame.
	void write_csv(const std::string &csv_name = "mcl_log.csv")
	{
		std::unordered_map<std::string, std::vector<double> > timings_ms;
		get_timings(timings_ms);
		std::lock_guard<std::mutex> guard(write_mutex);
		int num_frames = frame_data.size();
		if (num_frames == 0) { return; } // no computation (excluding init)

		std::stringstream ss;

		// Write header
		ss << "frame";
		std::unordered_map<std::string, std::vector<double> >::iterator it = timings_ms.begin();
		for (; it != timings_ms.end(); ++it)
		{
			if (it->first == "elapsed_ms") { ss << ",elapsed_s"; }
			else { ss << "," << it->first; }
		}

		ss << std::fixed << std::setprecision(6);

		// Loop frames (including init) and write timings
		// Special: if runtime_ms, store sum as seconds for precision purposes.
		double total_runtime_s = 0;
		for (int curr_write_frame = -1; curr_write_frame < num_frames; ++curr_write_frame)
		{
			ss << "\n" << curr_write_frame;
			for (it = timings_ms.begin(); it != timings_ms.end(); ++it)
			{
				double ms = (curr_write_frame==-1) ?
					init_data.timings_ms[it->first] :
					timings_ms[it->first][curr_write_frame];
				if (it->first == "elapsed_ms")
				{
					total_runtime_s += 0.001 * ms;
					ss << ',' << total_runtime_s;
				}
				else { ss << ',' << ms; }
			}
		}

		// Save to file
		std::ofstream ofout(csv_name.c_str());
		if (!ofout.good()) { return; }
		ofout << ss.str();
		ofout.close();
	}

protected:

	typedef std::unordered_map<std::string, MicroTimer>::iterator TimerIt;
	typedef std::unordered_map<std::string, double>::iterator DoubleIt;
	int curr_frame; // -1 = initialization
	std::mutex write_mutex;

	struct PerFrameData
	{
		// Starts a timer with the current label.
		void start_timer(const std::string &l)
		{
			TimerIt it = timers.emplace(l, MicroTimer()).first;
			it->second.reset(); // reset if the timer already exists
		}

		// Stops a timer with the current label.
		// If timings already exist for this label, they are summed.
		// Erases running timer when done
		void stop_timer(const std::string &l)
		{
			TimerIt it = timers.find(l);
			if (it == timers.end()) { return; } // timer never started
			DoubleIt it2 = timings_ms.emplace(l, 0).first;
			it2->second += it->second.elapsed_ms();
			timers.erase(it); // remove running timer
		}

		void stop_all_timers()
		{
			int num_timers = 0;
			std::vector<std::string> timer_names;
			for (TimerIt it = timers.begin(); it != timers.end(); ++it, ++num_timers) {
				timer_names.emplace_back(it->first);
			}
			for (int i=0; i<num_timers; ++i) { stop_timer(timer_names[i]); }
		}

		void add_time(const std::string &l, double sec)
		{
			DoubleIt it2 = timings_ms.emplace(l, 0).first;
			it2->second += sec;
		}

		std::unordered_map<std::string, double> timings_ms; // stopped timers
		std::unordered_map<std::string, MicroTimer> timers; // running timers
	};

	PerFrameData init_data;
	std::vector<PerFrameData> frame_data;
};

class FunctionTimerWrapper
{
public:
	const std::string f;
	static void reset_all() { mcl::Logger::get().clear(); }
	static void increment_frame() { mcl::Logger::get().increment_frame(); }
	static void stop_frame() { mcl::Logger::get().stop_frame(); }
	static void add_runtime(double sec, int iter) { mcl::Logger::get().add_runtime_s(sec, iter); }
	static void write_csv(const char *csv) { mcl::Logger::get().write_csv(std::string(csv)); }
	FunctionTimerWrapper(const char *f_) : f(f_) { mcl::Logger::get().start_timer(f); }
	~FunctionTimerWrapper() { mcl::Logger::get().stop_timer(f); }
};

} // ns mcl

#endif
