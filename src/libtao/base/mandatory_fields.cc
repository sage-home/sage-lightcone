#include "mandatory_fields.hh"
#include <algorithm>
#include <cctype>

namespace tao
{

// Case-insensitive string comparison
bool iequals(const std::string& a, const std::string& b)
{
    if (a.size() != b.size())
    {
        return false;
    }
    return std::equal(a.begin(), a.end(), b.begin(), [](char a, char b) {
        return std::tolower(static_cast<unsigned char>(a)) ==
               std::tolower(static_cast<unsigned char>(b));
    });
}

std::set<std::string> get_mandatory_sage_fields()
{
    return {sage_fields::POSX,
            sage_fields::POSY,
            sage_fields::POSZ,
            sage_fields::VELX,
            sage_fields::VELY,
            sage_fields::VELZ,
            sage_fields::SNAPNUM,
            sage_fields::GALAXY_INDEX,
            sage_fields::CENTRAL_GALAXY_INDEX,
            sage_fields::SAGE_TREE_INDEX,
            sage_fields::MERGE_TYPE,
            sage_fields::MERGE_INTO_ID,
            sage_fields::MERGE_INTO_SNAPNUM};
}

std::set<std::string> get_kdtree_required_fields()
{
    // Subset of mandatory fields that are absolutely critical for KD-tree
    // operation
    return {sage_fields::POSX, sage_fields::POSY, sage_fields::POSZ, sage_fields::SNAPNUM};
}

bool validate_sage_fields(const std::vector<std::string>& available_fields,
                          std::vector<std::string>& missing_fields)
{
    auto mandatory = get_mandatory_sage_fields();
    missing_fields.clear();

    for (const auto& required : mandatory)
    {
        bool found = false;
        for (const auto& available : available_fields)
        {
            if (iequals(required, available))
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            missing_fields.push_back(required);
        }
    }

    return missing_fields.empty();
}

} // namespace tao
