#ifndef tao_base_batch_hh
#define tao_base_batch_hh

#include <boost/any.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/map.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <libhpc/debug/assert.hh>
#include <libhpc/debug/except.hh>
#include <libhpc/system/assign.hh>
#include <libhpc/system/matrix.hh>
#include <libhpc/system/view.hh>
#include <string>
#include <unordered_map>

namespace tao {

template <class T> class batch {
public:
  enum field_value_type {
    STRING,
    DOUBLE,
    INTEGER,
    UNSIGNED,
    UNSIGNED_LONG,
    LONG_LONG,
    UNSIGNED_LONG_LONG,
    FLOAT,
    FIELD_VALUE_TERMINAL
  };

  enum field_rank_type { ATTRIBUTE, SCALAR, VECTOR };

  typedef T real_type;
  typedef std::tuple<boost::any, field_rank_type, field_value_type> field_type;

  typedef boost::mpl::map<
      boost::mpl::pair<std::string, boost::mpl::int_<STRING>>,
      boost::mpl::pair<double, boost::mpl::int_<DOUBLE>>,
      boost::mpl::pair<int, boost::mpl::int_<INTEGER>>,
      boost::mpl::pair<unsigned, boost::mpl::int_<UNSIGNED>>,
      boost::mpl::pair<unsigned long, boost::mpl::int_<UNSIGNED_LONG>>,
      boost::mpl::pair<long long, boost::mpl::int_<LONG_LONG>>,
      boost::mpl::pair<unsigned long long,
                       boost::mpl::int_<UNSIGNED_LONG_LONG>>,
      boost::mpl::pair<float, boost::mpl::int_<FLOAT>>>
      type_map;

public:
  batch() : _max_size(1000), _size(0), _mask(1000) {}

  batch(batch const &src)
      : _max_size(src._max_size), _size(src._size), _mask(src._mask) {
    for (const auto &item : src._fields) {
      const auto &name = item.first;
      const auto &src_field = item.second;
      auto &field = _fields[name];
      std::get<1>(field) = std::get<1>(src_field);
      std::get<2>(field) = std::get<2>(src_field);
      const auto &val = std::get<0>(src_field);
      if (std::get<1>(field) == ATTRIBUTE) {
        switch (std::get<2>(field)) {
        case STRING:
          std::get<0>(field) = boost::any_cast<std::string>(val);
          break;
        case DOUBLE:
          std::get<0>(field) = boost::any_cast<double>(val);
          break;
        case INTEGER:
          std::get<0>(field) = boost::any_cast<int>(val);
          break;
        case UNSIGNED:
          std::get<0>(field) = boost::any_cast<unsigned>(val);
          break;
        case UNSIGNED_LONG:
          std::get<0>(field) = boost::any_cast<unsigned long>(val);
          break;
        case LONG_LONG:
          std::get<0>(field) = boost::any_cast<long long>(val);
          break;
        case UNSIGNED_LONG_LONG:
          std::get<0>(field) = boost::any_cast<unsigned long long>(val);
          break;
        case FIELD_VALUE_TERMINAL:
          break;
        };
      } else if (std::get<1>(field) == SCALAR) {
        switch (std::get<2>(field)) {
        case STRING:
          std::get<0>(field) = new std::vector<std::string>(
              *boost::any_cast<std::vector<std::string> *>(val));
          break;
        case DOUBLE:
          std::get<0>(field) = new std::vector<double>(
              *boost::any_cast<std::vector<double> *>(val));
          break;
        case INTEGER:
          std::get<0>(field) =
              new std::vector<int>(*boost::any_cast<std::vector<int> *>(val));
          break;
        case UNSIGNED:
          std::get<0>(field) = new std::vector<unsigned>(
              *boost::any_cast<std::vector<unsigned> *>(val));
          break;
        case UNSIGNED_LONG:
          std::get<0>(field) = new std::vector<unsigned long>(
              *boost::any_cast<std::vector<unsigned long> *>(val));
          break;
        case LONG_LONG:
          std::get<0>(field) = new std::vector<long long>(
              *boost::any_cast<std::vector<long long> *>(val));
          break;
        case UNSIGNED_LONG_LONG:
          std::get<0>(field) = new std::vector<unsigned long long>(
              *boost::any_cast<std::vector<unsigned long long> *>(val));
          break;
        case FIELD_VALUE_TERMINAL:
          break;
        };
      } else {
        switch (std::get<2>(field)) {
        case STRING:
          std::get<0>(field) = new hpc::matrix<std::string>(
              *boost::any_cast<hpc::matrix<std::string> *>(val));
          break;
        case DOUBLE:
          std::get<0>(field) = new hpc::matrix<double>(
              *boost::any_cast<hpc::matrix<double> *>(val));
          break;
        case INTEGER:
          std::get<0>(field) =
              new hpc::matrix<int>(*boost::any_cast<hpc::matrix<int> *>(val));
          break;
        case UNSIGNED:
          std::get<0>(field) = new hpc::matrix<unsigned>(
              *boost::any_cast<hpc::matrix<unsigned> *>(val));
          break;
        case UNSIGNED_LONG:
          std::get<0>(field) = new hpc::matrix<unsigned long>(
              *boost::any_cast<hpc::matrix<unsigned long> *>(val));
          break;
        case LONG_LONG:
          std::get<0>(field) = new hpc::matrix<long long>(
              *boost::any_cast<hpc::matrix<long long> *>(val));
          break;
        case UNSIGNED_LONG_LONG:
          std::get<0>(field) = new hpc::matrix<unsigned long long>(
              *boost::any_cast<hpc::matrix<unsigned long long> *>(val));
          break;
        case FIELD_VALUE_TERMINAL:
          break;
        };
      }
    }
  }

  ~batch() { clear(); }

  std::unordered_map<std::string, field_type> getfields() const {
    return _fields;
  };

  batch &operator=(batch const &src) {
    _max_size = src._max_size;
    _size = src._size;
    hpc::assign(_mask, src._mask);

    for (const auto &item : src._fields) {
      const auto &name = item.first;
      const auto &src_field = item.second;
      auto &field = _fields[name];
      std::get<1>(field) = std::get<1>(src_field);
      std::get<2>(field) = std::get<2>(src_field);
      const auto &val = std::get<0>(src_field);
      if (std::get<1>(field) == ATTRIBUTE) {
        switch (std::get<2>(field)) {
        case STRING:
          std::get<0>(field) = boost::any_cast<std::string>(val);
          break;
        case DOUBLE:
          std::get<0>(field) = boost::any_cast<double>(val);
          break;
        case INTEGER:
          std::get<0>(field) = boost::any_cast<int>(val);
          break;
        case UNSIGNED:
          std::get<0>(field) = boost::any_cast<unsigned>(val);
          break;
        case UNSIGNED_LONG:
          std::get<0>(field) = boost::any_cast<unsigned long>(val);
          break;
        case LONG_LONG:
          std::get<0>(field) = boost::any_cast<long long>(val);
          break;
        case UNSIGNED_LONG_LONG:
          std::get<0>(field) = boost::any_cast<unsigned long long>(val);
          break;
        case FIELD_VALUE_TERMINAL:
          break;
        };
      } else if (std::get<1>(field) == SCALAR) {
        switch (std::get<2>(field)) {
        case STRING:
          std::get<0>(field) = new std::vector<std::string>(
              *boost::any_cast<std::vector<std::string> *>(val));
          break;
        case DOUBLE:
          std::get<0>(field) = new std::vector<double>(
              *boost::any_cast<std::vector<double> *>(val));
          break;
        case INTEGER:
          std::get<0>(field) =
              new std::vector<int>(*boost::any_cast<std::vector<int> *>(val));
          break;
        case UNSIGNED:
          std::get<0>(field) = new std::vector<unsigned>(
              *boost::any_cast<std::vector<unsigned> *>(val));
          break;
        case UNSIGNED_LONG:
          std::get<0>(field) = new std::vector<unsigned long>(
              *boost::any_cast<std::vector<unsigned long> *>(val));
          break;
        case LONG_LONG:
          std::get<0>(field) = new std::vector<long long>(
              *boost::any_cast<std::vector<long long> *>(val));
          break;
        case UNSIGNED_LONG_LONG:
          std::get<0>(field) = new std::vector<unsigned long long>(
              *boost::any_cast<std::vector<unsigned long long> *>(val));
          break;
        case FIELD_VALUE_TERMINAL:
          break;
        };
      } else {
        switch (std::get<2>(field)) {
        case STRING:
          std::get<0>(field) = new hpc::matrix<std::string>(
              *boost::any_cast<hpc::matrix<std::string> *>(val));
          break;
        case DOUBLE:
          std::get<0>(field) = new hpc::matrix<double>(
              *boost::any_cast<hpc::matrix<double> *>(val));
          break;
        case INTEGER:
          std::get<0>(field) =
              new hpc::matrix<int>(*boost::any_cast<hpc::matrix<int> *>(val));
          break;
        case UNSIGNED:
          std::get<0>(field) = new hpc::matrix<unsigned>(
              *boost::any_cast<hpc::matrix<unsigned> *>(val));
          break;
        case UNSIGNED_LONG:
          std::get<0>(field) = new hpc::matrix<unsigned long>(
              *boost::any_cast<hpc::matrix<unsigned long> *>(val));
          break;
        case LONG_LONG:
          std::get<0>(field) = new hpc::matrix<long long>(
              *boost::any_cast<hpc::matrix<long long> *>(val));
          break;
        case UNSIGNED_LONG_LONG:
          std::get<0>(field) = new hpc::matrix<unsigned long long>(
              *boost::any_cast<hpc::matrix<unsigned long long> *>(val));
          break;
        case FIELD_VALUE_TERMINAL:
          break;
        };
      }
    }

    return *this;
  }

  void clear() {
    _size = 0;
    for (auto &item : _fields) {
      auto &field = item.second;
      switch (std::get<1>(field)) {
      case SCALAR:
        _del_scalar(std::get<0>(field), std::get<2>(field));
        break;
      case VECTOR:
        _del_vector(std::get<0>(field), std::get<2>(field));
        break;
      default:
        break;
      };
    }
    _fields.clear();
    boost::fill(_mask, false);
  }

  void set_max_size(unsigned size) { _max_size = size; }

  void set_size(unsigned size) {
    if (size == 0 && _max_size == 0)
      _max_size = 1000;
    if (size > _max_size)
      _max_size = size;
    _size = size;
    _mask.resize(_size);
    boost::fill(_mask, false);
  }

  void update_size() {
    _size = std::numeric_limits<unsigned>::max();
    for (const auto &item : _fields) {
      const auto &field = item.second;
      const auto &val = std::get<0>(field);
      if (std::get<1>(field) == ATTRIBUTE)
        continue;
      switch (std::get<2>(field)) {
      case STRING:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size, boost::any_cast<std::vector<std::string> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size,
              boost::any_cast<hpc::matrix<std::string> *>(val)->n_rows());
        break;
      case DOUBLE:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size, boost::any_cast<std::vector<double> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size, boost::any_cast<hpc::matrix<double> *>(val)->n_rows());
        break;
      case INTEGER:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size, boost::any_cast<std::vector<int> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size, boost::any_cast<hpc::matrix<int> *>(val)->n_rows());
        break;
      case UNSIGNED:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size, boost::any_cast<std::vector<unsigned> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size, boost::any_cast<hpc::matrix<unsigned> *>(val)->n_rows());
        break;
      case UNSIGNED_LONG:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size,
              boost::any_cast<std::vector<unsigned long> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size,
              boost::any_cast<hpc::matrix<unsigned long> *>(val)->n_rows());
        break;
      case LONG_LONG:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size, boost::any_cast<std::vector<long long> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size, boost::any_cast<hpc::matrix<long long> *>(val)->n_rows());
        break;
      case UNSIGNED_LONG_LONG:
        if (std::get<1>(field) == SCALAR)
          _size = std::min<unsigned>(
              _size,
              boost::any_cast<std::vector<unsigned long long> *>(val)->size());
        else
          _size = std::min<unsigned>(
              _size, boost::any_cast<hpc::matrix<unsigned long long> *>(val)
                         ->n_rows());
        break;
      case FIELD_VALUE_TERMINAL:
        break;
      };
    }

    // Make sure we actually found something.
    if (_size == std::numeric_limits<unsigned>::max())
      _size = 0;

    _mask.resize(_size);
    boost::fill(_mask, false);
  }

  unsigned size() const { return _size; }

  unsigned max_size() const { return _max_size; }

  bool has_field(const std::string &name) const {
    return _fields.find(name) != _fields.end();
  }

  field_value_type get_field_type(const std::string &name) const {
    return std::get<2>(_fields.at(name));
  }

  template <class U>
  void set_attribute(const std::string &name, const U &value) {
    field_type &field = _fields[name];
    std::get<0>(field) = value;
    std::get<1>(field) = ATTRIBUTE;
    std::get<2>(field) =
        (field_value_type)boost::mpl::at<type_map, U>::type::value;
  }

  template <class U>
  typename hpc::view<std::vector<U>> set_scalar(const std::string &name) {
    ASSERT(_max_size, "Cannot set fields on batches with zero maximum size.");
    field_type &field = _fields[name];
    boost::any &val = std::get<0>(field);
    if (val.empty()) {
      val = new std::vector<U>(_max_size);
      std::get<1>(field) = (field_rank_type)SCALAR;
      std::get<2>(field) =
          (field_value_type)boost::mpl::at<type_map, U>::type::value;
    }
    return *boost::any_cast<std::vector<U> *>(val);
  }

  void set_scalar(const std::string &name, field_value_type type) {
    ASSERT(_max_size, "Cannot set fields on batches with zero maximum size.");
    field_type &field = _fields[name];
    boost::any &val = std::get<0>(field);
    if (val.empty()) {
      switch (type) {
      case STRING:
        val = new std::vector<std::string>(_max_size);
        break;
      case DOUBLE:
        val = new std::vector<double>(_max_size);
        break;
      case INTEGER:
        val = new std::vector<int>(_max_size);
        break;
      case UNSIGNED:
        val = new std::vector<unsigned>(_max_size);
        break;
      case UNSIGNED_LONG:
        val = new std::vector<unsigned long>(_max_size);
        break;
      case LONG_LONG:
        val = new std::vector<long long>(_max_size);
        break;
      case UNSIGNED_LONG_LONG:
        val = new std::vector<unsigned long long>(_max_size);
        break;
      case FIELD_VALUE_TERMINAL:
        break;
      default:
        ASSERT(0);
      }
      std::get<1>(field) = (field_rank_type)SCALAR;
      std::get<2>(field) = type;
    }
  }

  template <class U>
  hpc::matrix<U> &set_vector(const std::string &name, size_t size) {
    ASSERT(_max_size, "Cannot set fields on batches with zero maximum size.");
    field_type &field = _fields[name];
    boost::any &val = std::get<0>(field);
    if (val.empty()) {
      val = new hpc::matrix<U>(_max_size, size);
      std::get<1>(field) = (field_rank_type)VECTOR;
      std::get<2>(field) =
          (field_value_type)boost::mpl::at<type_map, U>::type::value;
    }
    return *boost::any_cast<hpc::matrix<U> *>(val);
  }

  template <class U> const U &attribute(const std::string &name) {
    return boost::any_cast<U &>(std::get<0>(field(name)));
  }

  template <class U>
  typename hpc::view<std::vector<U>> scalar(const std::string &name) {
    ASSERT(std::get<1>(field(name)) == SCALAR,
           "Assertion in batch object. Requesting "
           "field \"",
           name, "\" as SCALAR, however field is ",
           ((std::get<1>(field(name)) == 0) ? "ATTRIBUTE" : "VECTOR"), ".");
    return *boost::any_cast<std::vector<U> *>(std::get<0>(field(name)));
  }

  template <class U>
  const typename hpc::view<std::vector<U>>
  scalar(const std::string &name) const {
    ASSERT(std::get<1>(field(name)) == SCALAR,
           "Assertion in batch object. Requesting "
           "field \"",
           name, "\" as SCALAR, however field is ",
           ((std::get<1>(field(name)) == 0) ? "ATTRIBUTE" : "VECTOR"), ".");
    return *boost::any_cast<std::vector<U> *>(std::get<0>(field(name)));
  }

  template <class U> hpc::matrix<U> &vector(const std::string &name) {
    ASSERT(std::get<1>(field(name)) == VECTOR,
           "Assertion in batch object. Requesting "
           "field \"",
           name, "\" as VECTOR, however field is ",
           ((std::get<1>(field(name)) == 0) ? "ATTRIBUTE" : "SCALAR"), ".");
    return *boost::any_cast<hpc::matrix<U> *>(std::get<0>(field(name)));
  }

  template <class U>
  const hpc::matrix<U> &vector(const std::string &name) const {
    return *boost::any_cast<hpc::matrix<U> *>(std::get<0>(field(name)));
  }

  field_type &field(const std::string &name) {
    auto it = _fields.find(name);
    EXCEPT(it != _fields.end(), "Field not found on batch object: ", name);
    return it->second;
  }

  const field_type &field(const std::string &name) const {
    auto it = _fields.find(name);
    EXCEPT(it != _fields.end(), "Field not found on batch object: ", name);
    return it->second;
  }

  bool masked(unsigned idx) const { return _mask[idx]; }

  void mask(unsigned idx) { _mask[idx] = true; }

  void unmask(unsigned idx) { _mask[idx] = false; }

  friend std::ostream &operator<<(std::ostream &strm, const batch &obj) {
    return strm;
  }

protected:
  void _del_scalar(boost::any &val, field_value_type type) {
    switch (type) {
    case STRING:
      delete boost::any_cast<std::vector<std::string> *>(val);
      break;
    case DOUBLE:
      delete boost::any_cast<std::vector<double> *>(val);
      break;
    case INTEGER:
      delete boost::any_cast<std::vector<int> *>(val);
      break;
    case UNSIGNED:
      delete boost::any_cast<std::vector<unsigned> *>(val);
      break;
    case UNSIGNED_LONG:
      delete boost::any_cast<std::vector<unsigned long> *>(val);
      break;
    case LONG_LONG:
      delete boost::any_cast<std::vector<long long> *>(val);
      break;
    case UNSIGNED_LONG_LONG:
      delete boost::any_cast<std::vector<unsigned long long> *>(val);
      break;
    case FIELD_VALUE_TERMINAL:
      break;
    default:
      break;
    };
  }

  void _del_vector(boost::any &val, field_value_type type) {
    switch (type) {
    case STRING:
      delete boost::any_cast<hpc::matrix<std::string> *>(val);
      break;
    case DOUBLE:
      delete boost::any_cast<hpc::matrix<double> *>(val);
      break;
    case INTEGER:
      delete boost::any_cast<hpc::matrix<int> *>(val);
      break;
    case UNSIGNED:
      delete boost::any_cast<hpc::matrix<unsigned> *>(val);
      break;
    case UNSIGNED_LONG:
      delete boost::any_cast<hpc::matrix<unsigned long> *>(val);
      break;
    case LONG_LONG:
      delete boost::any_cast<hpc::matrix<long long> *>(val);
      break;
    case UNSIGNED_LONG_LONG:
      delete boost::any_cast<hpc::matrix<unsigned long long> *>(val);
      break;
    case FIELD_VALUE_TERMINAL:
      break;
    default:
      break;
    };
  }

protected:
  unsigned _max_size, _size;
  std::unordered_map<std::string, field_type> _fields;
  std::vector<bool> _mask;
};

} // namespace tao

#endif
