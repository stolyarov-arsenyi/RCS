#ifndef  RCS_STL_H
#define  RCS_STL_H


class STL
{
    using Face = std::array <float, 12>;

    using Mesh = std::vector <Face>;


    static auto parse_ascii (std::ifstream & file) -> Mesh
    {
        std::stringstream stream;

        stream << file.rdbuf();

        std::string string = stream.str();


        std::string rx_exponent = "e[+-]?[0-9]+";

        std::string rx_mantissa [] = { "[+-]?[0-9]+[.]?[0-9]*",
                                       "[+-]?[0-9]*[.]?[0-9]+" };

        std::string rx_float = "(?:" "(?:" + rx_mantissa[0] + ")" "|"
                                     "(?:" + rx_mantissa[1] + ")" ")" "(?:" + rx_exponent + ")?";

        std::string rx_space = "\\s+";

        std::string rx_point = "(" + rx_float + " " + rx_float + " " + rx_float + ")";

        std::string rx_facet = "facet normal point outer loop vertex point vertex point vertex point endloop endfacet";

        rx_facet = std::regex_replace(rx_facet, std::regex(R"(\s)"), rx_space);

        rx_facet = std::regex_replace(rx_facet, std::regex("point"), rx_point);

        Mesh mesh;

        std::smatch match;

        while (std::regex_search(string, match, std::regex(rx_facet)))
        {
            Face face {};

            std::stringstream(match[1]) >> face[0 + 0] >> face[0 + 1] >> face[0 + 2];
            std::stringstream(match[2]) >> face[3 + 0] >> face[3 + 1] >> face[3 + 2];
            std::stringstream(match[3]) >> face[6 + 0] >> face[6 + 1] >> face[6 + 2];
            std::stringstream(match[4]) >> face[9 + 0] >> face[9 + 1] >> face[9 + 2];

            mesh.push_back(face);

            string = match.suffix();
        }

        return mesh;
    }

    static auto parse_binary (std::ifstream & file) -> Mesh
    {
        file.seekg(80, std::ifstream::beg);


        std::uint32_t mesh_size;

        file.read((char *) & mesh_size, 4);


        Mesh mesh(mesh_size);

        for (auto & face : mesh)
        {
            file.read((char *) face.data(), 48);

            file.seekg(2, std::ifstream::cur);
        }

        return mesh;
    }

public:

    static auto parse (const std::string & name) -> Mesh
    {
        char head [5];

        std::ifstream file(name);

        if (file.read(head, 5))
        {
            if (std::string(head, 5) == "solid")

                return parse_ascii(file);

            return parse_binary(file);
        }

        std::cout << "error while reading file \"" << name << "\"" << std::endl;

        return {};
    }
};


#endif //RCS_STL_H
