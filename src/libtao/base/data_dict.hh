#ifndef tao_base_data_dict_hh
#define tao_base_data_dict_hh
#include <string>
#include "batch.hh"
#include "types.hh"
#include "xml_dict.hh"
#include "libhpc.hh"

/*
 * Dictionary of column names
 */

namespace tao {

	static const char* OBJECTTYPE = "objecttype";
    static const char* POSX = "posx";
    static const char* POSY = "posy";
    static const char* POSZ = "posz";
    static const char* VELX = "velx";
    static const char* VELY = "vely";
    static const char* VELZ = "velz";
    static const char* GLOBALINDEXID = "globalindexid";
    static const char* SNAPNUM = "snapnum";

    typedef union ubuftype {
        float f4;
        double f8;
        int32_t i4;
        int64_t i8;
    } ANY;

	class data_dict_field {
	public:
		data_dict_field()
		{

		}
		std::string _name;
		std::string _label;
        std::string _description;
		std::string _order;
		std::string _units;
		std::string _group;
        std::string _stats_minimum;
        std::string _stats_maximum;
        std::string _datatype;
        std::string _alias;
		typename tao::batch<real_type>::field_value_type _type;
		hpc::h5::datatype _mtype;
		hpc::h5::datatype _ftype;
        hsize_t _offset;
	};

	class data_dict {
	public:
	    data_dict() {

	    }
        static void setFieldStats(std::vector<tao::data_dict_field>& fields,std::string name, long long min, long long max);
        static void setFieldStats(std::vector<tao::data_dict_field>& fields,std::string name, int min, int max);
        static void setFieldStats(std::vector<tao::data_dict_field>& fields,std::string name, float min, float max);
        static xml_dict getSidecar(const std::string& fn, const char* rootnode="/settings");
	    static void saveSidecar(const std::string& fn, std::vector<tao::data_dict_field>& fields, bool as_settings=true);
        static std::vector<tao::data_dict_field> getFieldsANY(xml_dict& h5xml, const char* rootnode="/settings/sageinput/Field");
        static herr_t getFieldsH5( hid_t g_id, const char *name, const H5L_info_t *info, void *op_data);
        static int32_t findField(std::vector<tao::data_dict_field>& fields, std::string name);
        static int32_t findLabel(std::vector<tao::data_dict_field>& fields, std::string label);
        static hpc::h5::derive getComposite(std::vector<tao::data_dict_field>& fields,hpc::h5::datatype& mdt, hpc::h5::datatype& fdt);
        static hpc::h5::derive getCompositeANY(std::vector<tao::data_dict_field>& fields,hpc::h5::datatype& mdt, hpc::h5::datatype& fdt);
	};

}

#endif
