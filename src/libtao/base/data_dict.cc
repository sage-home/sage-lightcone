#include "data_dict.hh"
#include "libhpc.hh"

namespace tao
{

herr_t data_dict::getFieldsH5(hid_t g_id, const char* name, const H5L_info_t* info, void* op_data)
{
    std::vector<std::string>* fields = (std::vector<std::string>*)op_data;
    std::string name_lowercase = name;
    std::transform(name_lowercase.begin(), name_lowercase.end(), name_lowercase.begin(), ::tolower);
    fields->push_back(name);
    return 0;
}

int32_t data_dict::findField(std::vector<tao::data_dict_field>& fields, std::string name)
{
    // I want findField to be case insensitive
    std::string myname = name;
    hpc::to_lower(myname);
    for (uint32_t i = 0; i < fields.size(); i++)
    {
        std::string field_str = fields[i]._name;
        hpc::to_lower(field_str);
        if (field_str == myname)
            return i;
    }
    return -1;
}
int32_t data_dict::findLabel(std::vector<tao::data_dict_field>& fields, std::string label)
{
    // I want findField to be case insensitive
    std::string mylabel = label;
    hpc::to_lower(mylabel);
    for (uint32_t i = 0; i < fields.size(); i++)
    {
        std::string label_str = fields[i]._label;
        hpc::to_lower(label_str);
        if (label_str == mylabel)
            return i;
    }
    return -1;
}
hpc::h5::derive data_dict::getComposite(std::vector<tao::data_dict_field>& fields,
                                        hpc::h5::datatype& mdt, hpc::h5::datatype& fdt)
{
    // Build composite data type
    int size = 0;
    for (uint32_t i = 0; i < fields.size(); i++)
    {
        typename tao::batch<real_type>::field_value_type type = fields[i]._type;
        switch (type)
        {
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
            std::cout << "datatype not handled" << std::endl;
            break;
        }
    }
    hpc::h5::derive der(size);
    for (uint32_t i = 0; i < fields.size(); i++)
    {
        der.add(fields[i]._mtype, fields[i]._offset, fields[i]._ftype, fields[i]._name);
    }
    der.commit(mdt, fdt);
    return der;
}
hpc::h5::derive data_dict::getCompositeANY(std::vector<tao::data_dict_field>& fields,
                                           hpc::h5::datatype& mdt, hpc::h5::datatype& fdt)
{
    // Build composite data type
    int size = 0;
    for (uint32_t i = 0; i < fields.size(); i++)
    {
        typename tao::batch<real_type>::field_value_type type = fields[i]._type;
        switch (type)
        {
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
            std::cout << "datatype not handled" << std::endl;
            break;
        }
    }
    hpc::h5::derive der(size);
    for (uint32_t i = 0; i < fields.size(); i++)
    {
        der.add(fields[i]._mtype, fields[i]._offset, fields[i]._ftype, fields[i]._name);
    }
    der.commit(mdt, fdt);
    return der;
}
void data_dict::setFieldStats(std::vector<tao::data_dict_field>& fields, std::string name,
                              long long min, long long max)
{
    int index = findField(fields, name);
    if (index != -1)
    {
        fields[index]._stats_minimum = std::to_string(min);
        fields[index]._stats_maximum = std::to_string(max);
    }
}
void data_dict::setFieldStats(std::vector<tao::data_dict_field>& fields, std::string name, int min,
                              int max)
{
    int index = findField(fields, name);
    if (index != -1)
    {
        fields[index]._stats_minimum = std::to_string(min);
        fields[index]._stats_maximum = std::to_string(max);
    }
}
void data_dict::setFieldStats(std::vector<tao::data_dict_field>& fields, std::string name,
                              float min, float max)
{
    int index = findField(fields, name);
    if (index != -1)
    {
        fields[index]._stats_minimum = std::to_string(min);
        fields[index]._stats_maximum = std::to_string(max);
    }
}
} // namespace tao
