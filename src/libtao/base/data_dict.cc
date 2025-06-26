#include "data_dict.hh"
#include "libhpc.hh"

namespace tao {

    xml_dict data_dict::getSidecar(const std::string& fn, const char *rootnode) {
        std::string sidecar_filename = fn;
        size_t len=sidecar_filename.length();
        size_t pos=sidecar_filename.rfind('.');
        if (pos!=std::string::npos)
            sidecar_filename.erase(pos, len-pos);
        sidecar_filename.append(".xml");
        tao::xml_dict sidecar_xml;
        sidecar_xml.read(sidecar_filename, rootnode);
        return sidecar_xml;
    }

    void data_dict::saveSidecar(const std::string& fn, std::vector<tao::data_dict_field>& fields, bool as_settings)
    {
        std::string sidecar_filename = fn;
        size_t len=sidecar_filename.length();
        size_t pos=sidecar_filename.rfind('.');
        if (pos!=std::string::npos)
            sidecar_filename.erase(pos, len-pos);
        sidecar_filename.append(".xml");
        std::ofstream xmlfile;
        xmlfile.open( sidecar_filename, std::fstream::out | std::fstream::trunc );
        if (as_settings)
            xmlfile<<"<settings>\n";
        xmlfile<<"  <sageinput>\n";
        for (unsigned iii = 0; iii < fields.size(); iii++) {
            std::string typestr;
            typename tao::batch<real_type>::field_value_type type = fields[iii]._type;
            switch (type) {
                case tao::batch<real_type>::DOUBLE: {
                    typestr = "float";
                }
                    break;
                case tao::batch<real_type>::INTEGER: {
                    typestr = "int";
                }
                    break;
                case tao::batch<real_type>::LONG_LONG: {
                    typestr = "long long";
                }
                    break;
                default:
                    std::cout << "datatype not handled"<<std::endl;
                    break;
            }
            xmlfile<<"    <Field Type=\""<<typestr<<"\"\n";
            xmlfile<<"      label=\""<<fields[iii]._label<<"\"\n";
            xmlfile<<"      description=\""<<fields[iii]._description<<"\"\n";
            xmlfile<<"      order=\""<<fields[iii]._order<<"\"\n";
            xmlfile<<"      units=\""<<fields[iii]._units<<"\"\n";
            if (!fields[iii]._alias.empty()) {
                xmlfile<<"      alias=\""<<fields[iii]._alias<<"\"\n";
            }
            if (!fields[iii]._stats_minimum.empty()) {
                xmlfile<<"      stats_minimum=\""<<fields[iii]._stats_minimum<<"\"\n";
            }
            if (!fields[iii]._stats_maximum.empty()) {
                xmlfile<<"      stats_maximum=\""<<fields[iii]._stats_maximum<<"\"\n";
            }
            xmlfile<<"      group=\""<<fields[iii]._group<<"\">"<<fields[iii]._name<<"</Field>\n";
        }
        xmlfile<<"  </sageinput>\n";
        if (as_settings)
            xmlfile<<"</settings>\n";
        xmlfile.close();
    }

    herr_t data_dict::getFieldsH5( hid_t g_id, const char *name, const H5L_info_t *info, void *op_data) {
        std::vector<std::string> *fields = (std::vector<std::string>*)op_data;
        std::string name_lowercase=name;
        std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
        fields->push_back(name);
        return 0;
    }

    std::vector<tao::data_dict_field> data_dict::getFieldsANY(xml_dict& h5xml, const char* rootnode)
    {
        xpath_node_set fields = h5xml.get_nodes(rootnode);
        int order=0;
        hsize_t offset=0;
        std::vector<tao::data_dict_field> result;
        result.resize(fields.size());
        for (const xpath_node *it = fields.begin(); it != fields.end(); ++it) {
            xml_node cur = it->node();

            std::string element = cur.name();
            std::string field_str = cur.text().as_string();
            //hpc::to_lower(field_str);
            std::string label = cur.attribute("label").value();
            std::string description = cur.attribute("description").value();
            std::string strorder = cur.attribute("order").value();
            std::string units = cur.attribute("units").value();
            std::string group = cur.attribute("group").value();
            std::string type_str = cur.attribute("Type").value();
            std::string alias;
            std::string stats_minimum;
            std::string stats_maximum;
            if (cur.attribute("alias")) {
                alias = cur.attribute("alias").value();
            }
            if (cur.attribute("stats_minimum")) {
                stats_minimum = cur.attribute("stats_minimum").value();
                std::cout << "FOUND stats_minimum["<<field_str<<"]="<<stats_minimum<<std::endl;
            }
            if (cur.attribute("stats_maximum")) {
                stats_maximum = cur.attribute("stats_maximum").value();
            }
            typename tao::batch<real_type>::field_value_type type;
            tao::data_dict_field field;
            field._name = field_str;
            field._label = label;
            field._description = description;
            field._order = strorder;
            field._units = units;
            field._group = group;
            field._alias = alias;
            field._stats_minimum = stats_minimum;
            field._stats_maximum = stats_maximum;
            field._datatype = type_str;
            if (type_str == "int" || type_str == "INT" || type_str == "short" || type_str == "SMALLINT") {
                field._type = tao::batch<real_type>::INTEGER;
                field._mtype = hpc::h5::datatype::native_int;
                field._ftype = hpc::h5::datatype::std_i32be;
                field._offset = offset;
                offset += sizeof(ANY);
            } else if (type_str == "long long" || type_str == "BIGINT") {
                field._type = tao::batch<real_type>::LONG_LONG;
                field._mtype = hpc::h5::datatype::native_llong;
                field._ftype = hpc::h5::datatype::std_i64be;
                field._offset = offset;
                offset += sizeof(ANY);
            } else if (type_str == "float" || type_str == "FLOAT4" || type_str == "FLOAT8") {
                field._type = tao::batch<real_type>::DOUBLE;
                field._mtype = hpc::h5::datatype::native_float;
                field._ftype = hpc::h5::datatype::ieee_f32be;
                field._offset = offset;
                offset += sizeof(ANY);
            } else {
                EXCEPT(0, "Unknown field type for field '", field_str, "': ", type_str);
            }

            result[order] = field;
            order++;

        }
        return result;
    }
    int32_t data_dict::findField(std::vector<tao::data_dict_field>& fields,std::string name) {
        // I want findField to be case insensitive
        std::string myname=name;
        hpc::to_lower(myname);
        for (uint32_t i=0;i<fields.size();i++) {
            std::string field_str=fields[i]._name;
            hpc::to_lower(field_str);
            if (field_str == myname)
                return i;
        }
        return -1;
    }
    int32_t data_dict::findLabel(std::vector<tao::data_dict_field>& fields,std::string label) {
        // I want findField to be case insensitive
        std::string mylabel=label;
        hpc::to_lower(mylabel);
        for (uint32_t i=0;i<fields.size();i++) {
            std::string label_str=fields[i]._label;
            hpc::to_lower(label_str);
            if (label_str == mylabel)
                return i;
        }
        return -1;
    }
    hpc::h5::derive data_dict::getComposite(std::vector<tao::data_dict_field>& fields, hpc::h5::datatype& mdt, hpc::h5::datatype& fdt)
    {
        // Build composite data type
        int size=0;
        for (uint32_t i=0;i<fields.size();i++) {
            typename tao::batch<real_type>::field_value_type type = fields[i]._type;
            switch (type) {
                case tao::batch<real_type>::DOUBLE: {
                    size += sizeof(float);
                }
                    break;
                case tao::batch<real_type>::INTEGER: {
                    size += sizeof(int32_t);
                }
                    break;
                case tao::batch<real_type>::LONG_LONG: {
                    size += sizeof(long long);
                }
                    break;
                default:
                    std::cout << "datatype not handled"<<std::endl;
                    break;
            }
        }
        hpc::h5::derive der(size);
        for (uint32_t i=0;i<fields.size();i++) {
            der.add(fields[i]._mtype, fields[i]._offset, fields[i]._ftype, fields[i]._name);
        }
        der.commit(mdt, fdt);
        return der;
    }
    hpc::h5::derive data_dict::getCompositeANY(std::vector<tao::data_dict_field>& fields, hpc::h5::datatype& mdt, hpc::h5::datatype& fdt)
    {
        // Build composite data type
        int size=0;
        for (uint32_t i=0;i<fields.size();i++) {
            typename tao::batch<real_type>::field_value_type type = fields[i]._type;
            switch (type) {
                case tao::batch<real_type>::DOUBLE: {
                    size += sizeof(ANY);
                }
                    break;
                case tao::batch<real_type>::INTEGER: {
                    size += sizeof(ANY);
                }
                    break;
                case tao::batch<real_type>::LONG_LONG: {
                    size += sizeof(ANY);
                }
                    break;
                default:
                    std::cout << "datatype not handled"<<std::endl;
                    break;
            }
        }
        hpc::h5::derive der(size);
        for (uint32_t i=0;i<fields.size();i++) {
            der.add(fields[i]._mtype, fields[i]._offset, fields[i]._ftype, fields[i]._name);
        }
        der.commit(mdt, fdt);
        return der;
    }
    void data_dict::setFieldStats(std::vector<tao::data_dict_field>& fields,std::string name, long long min, long long max) {
        int index = findField(fields,name);
        if (index!=-1) {
            fields[index]._stats_minimum = std::to_string(min);
            fields[index]._stats_maximum = std::to_string(max);
        }
    }
    void data_dict::setFieldStats(std::vector<tao::data_dict_field>& fields,std::string name, int min, int max) {
        int index = findField(fields,name);
        if (index!=-1) {
            fields[index]._stats_minimum = std::to_string(min);
            fields[index]._stats_maximum = std::to_string(max);
        }
    }
    void data_dict::setFieldStats(std::vector<tao::data_dict_field>& fields,std::string name, float min, float max) {
        int index = findField(fields,name);
        if (index!=-1) {
            fields[index]._stats_minimum = std::to_string(min);
            fields[index]._stats_maximum = std::to_string(max);
        }
    }
}
