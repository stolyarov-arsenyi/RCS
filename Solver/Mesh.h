#ifndef RCS_MESH_MESH_H
#define RCS_MESH_MESH_H


template <class X>

struct Mesh
{
    std::vector <std::shared_ptr <Vector <Re <X>>>> vertices;

    std::vector <Edge <X>> edges;


    auto add_vertex (const Vector <Re <X>> & vertex) -> std::shared_ptr <Vector <Re <X>>>
    {
        for (const auto & v : vertices)

            if (* v == vertex)

                return v;

        return * vertices.insert(vertices.end(), std::make_shared <Vector <Re <X>>> (vertex));
    }

public:

    auto edge (std::size_t i) const -> const Edge <X> &
    {
        return edges[i];
    }

    auto size () const -> std::size_t
    {
        return edges.size();
    }


    static auto load_from_stl (const std::string & name) -> Mesh
    {
        struct Triangle
        {
            Vector <Re <X>> vertices [3];

            auto v (int i) const -> const Vector <Re <X>> &
            {
                return vertices[(3 + i % 3) % 3];
            }
        };

        auto data = STL::parse(name);

        Mesh mesh;

        for (std::size_t o = 0; o < data.size(); o ++)
        {
            Triangle outer = { { { data[o][3 + 0], data[o][3 + 1], data[o][3 + 2] },
                                 { data[o][3 + 3], data[o][3 + 4], data[o][3 + 5] },
                                 { data[o][3 + 6], data[o][3 + 7], data[o][3 + 8] } } };

            for (std::size_t i = o + 1; i < data.size(); i ++)
            {
                Triangle inner = { { { data[i][3 + 0], data[i][3 + 1], data[i][3 + 2] },
                                     { data[i][3 + 3], data[i][3 + 4], data[i][3 + 5] },
                                     { data[i][3 + 6], data[i][3 + 7], data[i][3 + 8] } } };

                for (int o = 0; o < 3; o ++)
                for (int i = 0; i < 3; i ++)
                {
                    if (outer.v(o + 1) == inner.v(i - 1)
                    &&  outer.v(o - 1) == inner.v(i + 1))
                    {
                        Face <X> face_outer (mesh.add_vertex(outer.v(o + 1)),
                                             mesh.add_vertex(outer.v(o - 1)),
                                             mesh.add_vertex(outer.v(o + 0)));

                        Face <X> face_inner (mesh.add_vertex(outer.v(o - 1)),
                                             mesh.add_vertex(outer.v(o + 1)),
                                             mesh.add_vertex(inner.v(i + 0)));

                        mesh.edges.emplace_back(face_outer, face_inner);
                    }
                }
            }
        }

        return mesh;
    }
};


#endif //RCS_MESH_MESH_H
