
#include "xml_dict.hh"
#include "libhpc/debug/except.hh"
#include <fstream>

using namespace pugi;

namespace tao {

xml_dict::xml_dict() : _sep(":"), _root(NULL) {}

xml_dict::xml_dict(xml_node root) : _sep(":"), _root(root) {}

xml_dict::xml_dict(std::istream &strm, const std::string &xpath_root)
    : _sep(":"), _root(0) {
  read(strm, xpath_root);
}

xml_dict::xml_dict(const xml_dict &src) : _sep(src._sep), _root(0) {
  _doc.reset(src._doc);
  if (src._root == src._doc)
    _root = _doc;
  else if (src._root)
    _root = _find_root(_doc, src._root.path());
}

xml_dict::~xml_dict() {}

void xml_dict::read(const std::string &filename, const std::string &path) {
  std::ifstream file(filename.c_str(), std::fstream::in);
  read(file, path, filename);
}

void xml_dict::read(std::istream &stream, const std::string &xpath_root,
                    const std::string &filename) {
  if (!_root) {
    xml_parse_result result = _doc.load(stream);
    EXCEPTAS(result, bad_xml, "Error: Couldn't read XML from file: ", filename);
    _root = _find_root(_doc, xpath_root);
  } else {
    _merge(stream, xpath_root, filename);
  }
}

bool xml_dict::has(const std::string &path) const {
  return _get_node(path, false);
}

xpath_node_set xml_dict::get_nodes(const std::string &xpath) const {
  return _root.select_nodes(xpath.c_str());
}

xml_node xml_dict::get_root() const { return _root; }

xml_node xml_dict::_find_root(xml_node &node,
                              const std::string &xpath_root) const {
  if (!xpath_root.empty()) {
    xpath_node root = node.select_single_node(xpath_root.c_str());
    ASSERT(root, "XPath root does not exist.");
    return root.node();
  } else
    return node;
}

void xml_dict::_merge(std::istream &stream, const std::string &path,
                      const std::string &filename) {
  xml_document doc;
  xml_parse_result result = doc.load(stream);
  EXCEPTAS(result, bad_xml, "Error: Couldn't read XML from file: ", filename);
  xml_node root = _find_root(doc, path);
  for (xml_node_iterator it = root.begin(); it != root.end(); ++it)
    _merge_node(_root, *it);
}

void xml_dict::_merge_node(xml_node merge_into, xml_node merge_from) {
  if (merge_from.type() == node_element) {
    xml_node child = merge_into.child(merge_from.name());
    if (child) {
      for (xml_node_iterator it = merge_from.begin(); it != merge_from.end();
           ++it)
        _merge_node(child, *it);
    } else
      merge_into.append_copy(merge_from);
  } else if (merge_from.type() == node_pcdata) {
    // Search for the first pcdata element.
    for (xml_node_iterator it = merge_into.begin(); it != merge_into.end();
         ++it) {
      if (it->type() == node_pcdata) {
        it->set_value(merge_from.value());
        return;
      }
    }
  }
}

xml_node xml_dict::_get_node(const std::string &path, bool except) const {
  std::string xpath = _xform_path(path);
  xpath_node node = _root.select_single_node(xpath.c_str());
  if (!node) {
    EXCEPTAS(!except, bad_option,
             "Error: Couldn't find option in XML dictionary: ", path);
    return (xml_node)0;
  }
  return node.node();
}

std::string xml_dict::_xform_path(const std::string &path) const {
  std::string new_path = std::string("./") + path;
  for (unsigned ii = 0; ii < new_path.size(); ++ii) {
    if (new_path[ii] == ':')
      new_path[ii] = '/';
  }
  return new_path;
}

template <> std::string xml_dict::_coerce(const std::string &value) const {
  return value;
}

} // namespace tao
