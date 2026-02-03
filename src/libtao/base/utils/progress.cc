#include "progress.hh"

namespace tao {

progress::progress() : _comm(hpc::mpi::comm::world) {}

void progress::start(double total, double running) {
  _done = false;
  _cache.resize(_comm.size());
  _history.resize(max_history);
  _history_pos = 0;
  _last_time.resize(_comm.size());
  std::fill(_last_time.begin(), _last_time.end(), 0.0);
  _total = total;
  _running = running;
  _all_time = 0.0;
  _start = std::chrono::high_resolution_clock::now();
  _sendt = std::thread(&progress::_sender, this);
  if (_comm.rank() == 0) {
    _recvt = std::thread(&progress::_receiver, this);
    _rept = std::thread(&progress::_reporter, this);
  }
}

void progress::stop() {
  value_type val = std::make_pair(-1.0, -1.0);
  _comm.send(val, 0, progress_tag);
  _done = true;
  if (_sendt.joinable())
    _sendt.join();
  if (_recvt.joinable())
    _recvt.join();
  if (_rept.joinable())
    _rept.join();
}

void progress::append(double val) {
  _lock.lock();
  _out.emplace_back(_get_time(), val);
  _lock.unlock();
}

void progress::_receiver() {
  int n_finished = 0;
  while (n_finished < _comm.size()) {
    MPI_Status stat;
    if (_comm.iprobe(stat, MPI_ANY_SOURCE, progress_tag)) {
      value_type val;
      _comm.recv((double *)&val, hpc::mpi::datatype::double_floating,
                 stat.MPI_SOURCE, 2, stat.MPI_TAG);
      _lock.lock();
      if (val.first < 0.0) {
        _last_time[stat.MPI_SOURCE] = -1.0;
        ++n_finished;
      } else {
        _cache[stat.MPI_SOURCE].push_back(val);
        _current += val.second;
      }
      _lock.unlock();
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }
  }
}

void progress::_sender() {
  while (!_done) {
    while (_out.size()) {
      value_type val;
      _lock.lock();
      val = _out.front();
      _out.pop_front();
      _lock.unlock();
      hpc::mpi::request req;
      _comm.isend(val, 0, req, progress_tag);
      while (!_done && !req.test()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
      }
      if (_done)
        break;
    }
  }
}

void progress::_reporter() {
  while (!_done) {
    std::this_thread::sleep_for(std::chrono::milliseconds(report_sleep_ms));

    // Lock here to avoid conflicts.
    _lock.lock();

    // Recalculate historical values.
    _update_history();

    // Update the running values, the current round values, and the remaining.
    _running += _current;
    _current = 0;
    double remaining = _total - _running;
    double perc = 100.0 * _running / _total;

    // Allow access to the data cache again.
    _lock.unlock();

    // Calculate a global rate and time-to-go.
    _all_time += (report_sleep_ms / 1000.0);
    double all_rate = _running / _all_time;
    double all_ttg = remaining / all_rate;

    // Take the average of the historical values. It's easy, and seems to
    // provide a pretty good estimate.
    double avg_history_rate = 0.0;
    unsigned cnt = 0;
    for (auto const &x : _history) {
      if (x.first > 0.0) {
        avg_history_rate += x.second;
        ++cnt;
      }
    }
    avg_history_rate /= (double)cnt;
    double avg_history_ttg = remaining / avg_history_rate;

    // Take some of the global average and some of the recent history average.
    double best_ttg = 0.5 * all_ttg + 0.5 * avg_history_ttg;

    // Display some output.
    LOG_PUSH_TAG("progress");
    LOGILN(_get_time(), ",", perc, "%,", avg_history_ttg, "s");
    LOG_POP_TAG("progress");
  }
}

void progress::_update_history() {
  // The way this routine operates is based on maintaining a list of received
  // values for each rank. When new values come in, we stash them in the source
  // rank's list and then call this routine. Each entry in the cache corresponds
  // to a 2D point, the coordinates of which are number of objects processed,
  // and when the measurement was taken.
  //
  // Once in here, we try to find as many points in time where we have an entry
  // in the cache from each rank available. When we find one such point we can
  // calculate the global instantaneous average rate of object processing, and
  // push that to our rate list. If we should ever find that any one of our
  // lists is exhausted, we can terminate this routine and wait for more entries
  // to arrive.
  //
  // `while( 1 )` is used as we don't know how many complete entries in the
  // cache we'll find, so it's easier to drop out once we exhaust one rank's
  // list.
  while (1) {
    value_type val;
    double lowest_time = std::numeric_limits<double>::max(), sum = 0.0;
    unsigned lowest_rank;
    unsigned n_finished = 0;
    for (unsigned ii = 0; ii < _cache.size(); ++ii) {
      if (_last_time[ii] < 0.0) {
        ++n_finished;
        continue;
      }

      // This is where we terminate, as we've hit the end of one of our cache
      // lists.
      if (!_cache[ii].size())
        return;

      val = _cache[ii].front();
      if (val.first < lowest_time) {
        lowest_time = val.first;
        lowest_rank = ii;
      }
      double duration = val.first - _last_time[ii];
      double rate = val.second / duration;
      sum += rate;
    }

    // We need this here to prevent deadlock; because we eliminate completed
    // ranks from the calculation the previous `return` will get skipped once
    // all ranks have completed.
    if (n_finished == static_cast<unsigned>(_comm.size()))
      return;

    val = _cache[lowest_rank].front();
    _last_time[lowest_rank] = val.first;
    _cache[lowest_rank].pop_front();
    _history[_history_pos++] = std::make_pair(lowest_time, sum);
    _history_pos = _history_pos % _history.size();
  }
}

double progress::_get_time() {
  auto cur = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = cur - _start;
  return diff.count();
}

} // namespace tao
