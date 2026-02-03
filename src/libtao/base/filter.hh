#ifndef tao_base_filter_hh
#define tao_base_filter_hh

#include "batch.hh"
#include "types.hh"
#include <boost/optional.hpp>
#include <string>

namespace tao {
using namespace hpc;

class filter_iterator;

/// Defines a filter over a single field of a data source. Uses a
/// minimum and maximum value to skip batch entries.
///
class filter {
public:
  typedef filter_iterator iterator;

public:
  /// Constcut an empty instance.
  ///
  filter();

  /// Test a single index within a batch object to see if it passes
  /// the filter's test.
  ///
  /// @param[in] bat The batch object.
  /// @param[in] idx The index of the entry within the batch.
  /// @returns bool True if the batch entry passes the filter.
  ///
  bool operator()(batch<real_type> const &bat, bool filter_field_enabled,
                  unsigned idx) const {
    if (bat.masked(idx))
      return false;
    if (!filter_field_enabled)
      return true;
    switch (_type) {
    case batch<real_type>::DOUBLE: {
      double val = bat.scalar<double>(_field_name)[idx];
      if (_min && boost::any_cast<double>(*_min) > val)
        return false;
      if (_max && boost::any_cast<double>(*_max) < val)
        return false;
    } break;

    case batch<real_type>::INTEGER: {
      int val = bat.scalar<int>(_field_name)[idx];
      if (_min && boost::any_cast<int>(*_min) > val)
        return false;
      if (_max && boost::any_cast<int>(*_max) < val)
        return false;
    } break;

    case batch<real_type>::UNSIGNED_LONG: {
      unsigned long val = bat.scalar<unsigned long>(_field_name)[idx];
      if (_min && boost::any_cast<unsigned long>(*_min) > val)
        return false;
      if (_max && boost::any_cast<unsigned long>(*_max) < val)
        return false;
    } break;

    case batch<real_type>::LONG_LONG: {
      long long val = bat.scalar<long long>(_field_name)[idx];
      if (_min && boost::any_cast<long long>(*_min) > val)
        return false;
      if (_max && boost::any_cast<long long>(*_max) < val)
        return false;
    } break;

    case batch<real_type>::UNSIGNED_LONG_LONG: {
      unsigned long long val = bat.scalar<unsigned long long>(_field_name)[idx];
      if (_min && boost::any_cast<unsigned long long>(*_min) > val)
        return false;
      if (_max && boost::any_cast<unsigned long long>(*_max) < val)
        return false;
    }

    default:
      break;
    }
    return true;
  }

  /// Set the field and its limits. The field name is given as a string, and
  /// must match a field stored on the batch object to be filtered. The minimum
  /// and maximum values are also given as string representations of their
  /// values.
  ///
  /// @param[in] name The field name.
  /// @param[in] min The minimum value.
  /// @param[in] max The maximum value.
  ///
  void set_field(const std::string &name, const std::string &min,
                 const std::string &max);

  /// Set the type of the field. This indicates how the filter should interpret
  /// the string values given in `set_field`.
  ///
  /// @param[in] type The field type.
  ///
  void set_type(batch<real_type>::field_value_type type);

  template <class T> void set_minimum(const T &val) {
    switch (_type) {
    case batch<real_type>::DOUBLE:
      _min = (double)val;
      break;

    case batch<real_type>::INTEGER:
      _min = (int)val;
      break;

    case batch<real_type>::UNSIGNED_LONG:
      _min = (unsigned long)val;
      break;

    case batch<real_type>::LONG_LONG:
      _min = (long long)val;
      break;

    case batch<real_type>::UNSIGNED_LONG_LONG:
      _min = (unsigned long long)val;
      break;

    default:
      break;
    }
  }

  template <class T> void set_maximum(const T &val) {
    switch (_type) {
    case batch<real_type>::DOUBLE:
      _max = (double)val;
      break;

    case batch<real_type>::INTEGER:
      _max = (int)val;
      break;

    case batch<real_type>::UNSIGNED_LONG:
      _max = (unsigned long)val;
      break;

    case batch<real_type>::LONG_LONG:
      _max = (long long)val;
      break;

    case batch<real_type>::UNSIGNED_LONG_LONG:
      _max = (unsigned long long)val;
      break;

    default:
      break;
    }
  }

  const std::string &field_name() const;

  template <class T> boost::optional<T> minimum() const {
    if (_min) {
      T val;
      switch (_type) {
      case batch<real_type>::DOUBLE:
        val = boost::any_cast<double>(*_min);
        break;

      case batch<real_type>::INTEGER:
        val = boost::any_cast<int>(*_min);
        break;

      case batch<real_type>::UNSIGNED_LONG:
        val = boost::any_cast<unsigned long>(*_min);
        break;

      case batch<real_type>::LONG_LONG:
        val = boost::any_cast<long long>(*_min);
        break;

      case batch<real_type>::UNSIGNED_LONG_LONG:
        val = boost::any_cast<unsigned long long>(*_min);
        break;
      default:
        return boost::optional<T>(0);
      }
      return boost::optional<T>(val);
    } else
      return boost::none;
  }

  template <class T> boost::optional<T> maximum() const {
    if (_max) {
      T val;
      switch (_type) {
      case batch<real_type>::DOUBLE:
        val = boost::any_cast<double>(*_max);
        break;

      case batch<real_type>::INTEGER:
        val = boost::any_cast<int>(*_max);
        break;

      case batch<real_type>::UNSIGNED_LONG:
        val = boost::any_cast<unsigned long>(*_max);
        break;

      case batch<real_type>::LONG_LONG:
        val = boost::any_cast<long long>(*_max);
        break;

      case batch<real_type>::UNSIGNED_LONG_LONG:
        val = boost::any_cast<unsigned long long>(*_max);
        break;

      default:
        return boost::optional<T>(0);
      }
      return boost::optional<T>(val);
    } else
      return boost::none;
  }

  /// Iterate over the entries of a batch object, skipping those
  /// entries that fail the filter test.
  ///
  /// @param[in] bat The batch object to be filtered.
  /// @returns An iterator at the beginning of the batch object.
  ///
  iterator begin(const tao::batch<real_type> &bat,
                 bool filter_field_enabled) const;

  /// Iterate over the entries of a batch object, skipping those
  /// entries that fail the filter test.
  ///
  /// @param[in] bat The batch object to be filtered.
  /// @returns An iterator at the end of the batch object.
  ///
  iterator end(const tao::batch<real_type> &bat,
               bool filter_field_enabled) const;

  void settletomask(tao::batch<real_type> &bat,
                    bool filter_field_enabled) const {
    unsigned int ii = 0;
    while (ii < bat.size()) {
      if (!(*this)(bat, filter_field_enabled, ii))
        bat.mask(ii);
      ++ii;
    }
  }
  bool calculated_fields_are_ready(bool yesno) const {
    bool fields_are_ready = true;
    if (!yesno) {
      if (ends_with(_field_name, "_apparent")) {
        fields_are_ready = false;
      }
      if (ends_with(_field_name, "_absolute")) {
        fields_are_ready = false;
      }
    }
    return fields_are_ready;
  }

  inline bool ends_with(std::string const &value,
                        std::string const &ending) const {
    if (ending.size() > value.size())
      return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  }

protected:
  std::string _field_name;
  std::string _min_str;
  std::string _max_str;
  boost::optional<boost::any> _min;
  boost::optional<boost::any> _max;
  batch<real_type>::field_value_type _type;
};

class filter_iterator
    : public boost::iterator_facade<filter_iterator, unsigned,
                                    std::forward_iterator_tag, unsigned> {
  friend class boost::iterator_core_access;

public:
  filter_iterator();

  filter_iterator(const filter *filt, const batch<real_type> *bat,
                  bool filter_field_enabled, bool done);

protected:
  void increment();

  bool equal(const filter_iterator &op) const;

  unsigned dereference() const;

  void _settle(bool field_filter_enabled);

protected:
  const bool _filter_field_enabled;
  const tao::filter *_filt;
  const tao::batch<real_type> *_bat;
  unsigned _idx;
};

} // namespace tao

#endif
