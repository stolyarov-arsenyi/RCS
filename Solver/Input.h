#ifndef RCS_SOLVER_INPUT_H
#define RCS_SOLVER_INPUT_H


struct Input : std::map <std::string, std::string>
{
    static auto parse (const std::string & name) -> Input
    {
        Input input;

        std::ifstream file(name);

        std::string string;

        while (std::getline(file, string))
        {
            std::stringstream line(string);

            std::string name;
            std::string tail;

            line >> name;

            if (name.empty())

                continue;

            std::getline(line, tail);

            input[name] = tail;
        }

        return input;
    }
};


#endif //RCS_SOLVER_INPUT_H
