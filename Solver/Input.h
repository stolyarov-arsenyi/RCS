#ifndef RCS_SOLVER_INPUT_H
#define RCS_SOLVER_INPUT_H


template <class X>

struct Input
{
    static auto parse (const std::string & name) -> std::map <std::string, std::vector <X>>
    {
        std::map <std::string, std::vector <X>> input;

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

                for (X val = range[0]; val < range[1]; val += range[2])

                    input[token].push_back(val);
        }

        if (input["wavelength"  ].empty()) input["wavelength"  ].push_back(1);
        if (input["polarization"].empty()) input["polarization"].push_back(0);
        if (input["azimuth"     ].empty()) input["azimuth"     ].push_back(0);
        if (input["altitude"    ].empty()) input["altitude"    ].push_back(0);
        if (input["threads"     ].empty()) input["threads"     ].push_back(1);
        if (input["block_size"  ].empty()) input["block_size"  ].push_back(1);
        if (input["save_matrix" ].empty()) input["save_matrix" ].push_back(0);

        return input;
    }
};


#endif //RCS_SOLVER_INPUT_H
