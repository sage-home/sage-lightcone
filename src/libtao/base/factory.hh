#ifndef tao_base_factory_hh
#define tao_base_factory_hh

#include "module.hh"
#include <boost/iterator/transform_iterator.hpp>
#include <libhpc/system/has.hh>
#include <list>
#include <map>
#include <string>

namespace tao {

template <class Backend> class factory {
public:
  typedef Backend backend_type;
  typedef module<backend_type> module_type;
  typedef module_type *(*factory_create_type)(const std::string &name,
                                              pugi::xml_node base);
  typedef typename std::list<module_type *>::iterator iterator;

public:
  factory() {}

  ~factory() {
    for (auto mod : _mods)
      delete mod;
  }

  void register_module(const std::string &name, factory_create_type create) {
    EXCEPT(!hpc::has(_facs, name),
           "Cannot add new science module factory, "
           "ID already exists: ",
           name);
    _facs.emplace(name, create);
  }

  module_type &create_module(const std::string &name,
                             const std::string &inst_name = std::string(),
                             pugi::xml_node base = (pugi::xml_node)0) {
    std::string _in;
    if (inst_name.empty())
      _in = name;
    else
      _in = inst_name;
#ifndef NEXCEPT
//        for( auto mod : _mods )
//            EXCEPT( mod->name() != _in, "Science module with name already
//            exists: ", _in );
#endif
    EXCEPT(hpc::has(_facs, name), "No science module exists with ID: ", name);
    module_type *mod = _facs.at(name)(_in, base);
    _mods.push_back(mod);
    return *mod;
  }

  module_type &
  create_module_without_xml(const std::string &name,
                            const std::string &inst_name = std::string(),
                            pugi::xml_node base = (pugi::xml_node)0) {
    std::string _in;
    if (inst_name.empty())
      _in = name;
    else
      _in = inst_name;
#ifndef NEXCEPT
//        for( auto mod : _mods )
//            EXCEPT( mod->name() != _in, "Science module with name already
//            exists: ", _in );
#endif
    EXCEPT(hpc::has(_facs, name), "No science module exists with ID: ", name);
    module_type *mod = _facs.at(name)(_in, base);
    _mods.push_back(mod);
    LOGDLN("Created module ", name, " with name ", _in);
    return *mod;
  }

  iterator begin() { return _mods.begin(); }

  iterator end() { return _mods.end(); }

  module_type *operator[](const std::string &name) {
    for (auto mod : _mods) {
      if (mod->name() == name)
        return mod;
    }
    EXCEPT(0, "Failed to locate module with name: ", name);
  }

protected:
  std::map<std::string, factory_create_type> _facs;
  std::list<module_type *> _mods;
};

} // namespace tao

#endif
