#ifndef RCS_SOLVER_INPUT_H
#define RCS_SOLVER_INPUT_H


template <class X>

struct Input : std::map <std::string, std::vector <X>>
{
    static auto parse (const std::string & name) -> Input
    {
        Input input;

        std::ifstream file(name);

        std::string string;

        while (std::getline(file, string))
        {
            std::stringstream line(string);

            std::string token;

            line >> token;

            std::vector <X> range;

            X value;

            while (line >> value)

                range.push_back(value);

            if (range.size() == 1)

                input[token].push_back(range[0]);

            if (range.size() == 3)
            {
                for (X val = range[0]; val < range[1]; val += range[2])

                    input[token].push_back(val);

                input[token].push_back(range[1]);
            }
        }

        if (input["threads"   ].empty()) input["threads"   ].push_back(1);
        if (input["block_size"].empty()) input["block_size"].push_back(1);

        if (input["wavelength"].empty()) input["wavelength"].push_back(1);

        if (input["pol_inc"].empty()) input["pol_inc"].push_back(0);
        if (input["azi_inc"].empty()) input["azi_inc"].push_back(0);
        if (input["alt_inc"].empty()) input["alt_inc"].push_back(0);

        if (input["pol_sca"].empty()) input["pol_sca"].push_back(0);
        if (input["azi_sca"].empty()) input["azi_sca"].push_back(0);
        if (input["alt_sca"].empty()) input["alt_sca"].push_back(0);

        return input;
    }
};


#endif //RCS_SOLVER_INPUT_H
