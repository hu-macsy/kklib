#include <kklib/metascheme.hpp>

std::vector<std::vector<std::vector<bool>>> read_metapath_schemes(const char* path)
{
    FILE* f = fopen(path, "r");
    assert(f != NULL);
    int scheme_num, state_num;
    assert(2 == fscanf(f, "%d %d", &scheme_num, &state_num));
    std::vector<std::vector<std::vector<bool>>> schemes(scheme_num);
    for (int sc_i = 0; sc_i < scheme_num; sc_i++)
    {
        int scheme_length;
        assert(1 == fscanf(f, "%d", &scheme_length));
        schemes[sc_i].resize(scheme_length);
        for (int l_i = 0; l_i < scheme_length; l_i++)
        {
            schemes[sc_i][l_i].resize(state_num);
            for (int st_i = 0; st_i < state_num; st_i++)
            {
                int st;
                assert(1 == fscanf(f, "%d", &st));
                schemes[sc_i][l_i][st_i] = (st != 0);
            }
        }
    }
    return schemes;
    fclose(f);
};

void write_metapath_schemes(std::vector<std::vector<std::vector<bool>>> schemes, const char* path)
{
    FILE* f = fopen(path, "w");
    assert(f != NULL);
    int scheme_num = schemes.size();
    int state_num = schemes[0][0].size();
    fprintf(f, "%d %d\n", scheme_num, state_num);
    for (int sc_i = 0; sc_i < scheme_num; sc_i++)
    {
        fprintf(f, "%d\n", (int)schemes[sc_i].size());
        for (auto trans : schemes[sc_i])
        {
            for (auto s : trans)
            {
                fprintf(f, "%d ", (int)s);
            }
            fprintf(f, "\n");
        }
    }
    fclose(f);
}

std::vector<std::vector<scheme_mask_t>> get_scheme_mask(std::vector<std::vector<std::vector<bool>>>& schemes)
{
    std::vector<std::vector<scheme_mask_t>> ret(schemes.size());

    for (size_t i = 0; i < schemes.size(); i++)
    {
        ret[i].resize(schemes[i].size());
        for (size_t j = 0; j < schemes[i].size(); j++)
        {
            scheme_mask_t val = 0;
            for (size_t k = 0; k < schemes[i][j].size(); k++)
            {
                if (schemes[i][j][k])
                {
                    val |= (1 << k);
                }
            }
            ret[i][j] = val;
        }
    }
    return ret;
}