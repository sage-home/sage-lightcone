#ifndef tao_base_utils_progress_hh
#define tao_base_utils_progress_hh

#include <boost/optional.hpp>
#include <chrono>
#include <libhpc/mpi/comm.hh>
#include <mutex>
#include <thread>

namespace tao
{

/// Distributed progress tracking utility.
///
/// Fires up threads to handle sending and receiving data in the background.
/// The algorithm uses a combination of a global average and an average of the
/// the most recent n calculations. Steps are taken to ensure we only use
/// values that have meaning on all ranks.
///
/// The rate is calculated by ranks each sending in how many objects they've
/// processed (objects in this case are galaxies). We keep an estimate of the
/// total number of objects to process in order to provide an estimated time
/// to completion.
///
class progress
{
public:
    int const progress_tag = 434359;   // chosen at random
    int const sleep_ms = 500;          // time to wait between recv checks
    int const report_sleep_ms = 10000; // time between outputs
    int const max_history = 10;        // how many recent history values to keep
    typedef std::pair<double, double> value_type;

public:
    progress();

    /// Begin tracking progress. This fires up the threads and clears all
    /// initial values.
    ///
    /// @param total The total number of objects expected.
    ///
    void start(double total, double running = 0);

    /// Terminate the progress tracking. This cleans up the threads.
    ///
    void stop();

    /// Append a value to the outgoing queue.
    ///
    /// @param val The value to send (number of objects).
    ///
    void append(double val);

protected:
    /// The receiver thread. This is launched during `start`.
    ///
    void _receiver();

    /// The sender thread. This is launched during `start`.
    ///
    void _sender();

    /// The reporter thread. This is launched during `start` and handles
    /// calculating the global average and displaying the remaining time.
    ///
    void _reporter();

    /// Recalculate the recent history values. Behaves a bit like
    /// a ciruclar buffer.
    ///
    void _update_history();

    /// Get the time elapsed since timing began.
    ///
    double _get_time();

protected:
    std::list<value_type> _out;
    std::vector<std::list<value_type>> _cache;
    std::vector<value_type> _history;
    unsigned _history_pos;
    std::vector<double> _last_time;
    double _current, _running, _total, _all_time;
    std::chrono::high_resolution_clock::time_point _start;
    bool _done;
    std::thread _recvt, _sendt, _rept;
    std::mutex _lock;
    hpc::mpi::comm _comm;
};

} // namespace tao

#endif
