#include <print>

#include "input_gmsh.h"
#include "output_silo.h"
#include "geom_mesh.h"
#include "graph.h"

int main(int argc, const char **argv)
{
    if (argc < 2) {
        std::println("Specify mesh file name");
        return 1;
    }

    /**********************************************
     * Mesh loading
     */
    /* Declare the mesh object */
    frico::mesh msh;

    /* Load from the GMSH .geo file */
    auto meshok = frico::load_from_gmsh(argv[1], msh);
    if (not meshok) {
        std::println(stderr, "Can't load {}", argv[1]);
        return 1;
    }

    /**********************************************
     * Mesh connectivity:
     */
    /* The connectivity between edges and triangles is
     * contained in msh.edge_neighbours:
     */
    for (size_t iedg = 0; iedg < msh.edge_neighbours.size(); iedg++) {
        /* iedg is not only the counter but also the edge number */
        /* T- is always the triangle with the smallest index, while
         * T+ is always the triangle with the biggest index. The index
         * of T- is a plain integer contained in .itminus, while .itplus
         * is a std::optional: If the edge is a boundary edge, .itplus
         * will be an empty optional. */
        const auto& en = msh.edge_neighbours[iedg]; // conn info in 'en'
        if (not en.itplus) {
            // If en.itplus is empty this
            // is a boundary edge: handle as needed
            continue;
        }
    
        assert(en.itminus < msh.triangles.size());
        auto Tminus_num = en.itminus;               // T- number
        auto Tminus = msh.triangles[Tminus_num];    // T+ triangle
        assert(en.itplus.value() < msh.triangles.size());
        auto Tplus_num = en.itplus.value();         // T+ number
        auto Tplus = msh.triangles[Tplus_num];      // T+ triangle
    }

    /* The connectivity between edges and vertices is
     * contained directly into the edge object */
    for (size_t iedg = 0; iedg < msh.edges.size(); iedg++) {
        /* iedg is not only the counter but also the edge number */
        const auto& edge = msh.edges[iedg];
        auto evtx0 = edge.iv0; /* Index of first vertex */
        auto evtx1 = edge.iv1; /* Index of second vertex */
    }

    /* Connectivity between triangles and vertices is
     * contained in the triangle objects and it is obtained
     * exactly as for edges
     */
    for (size_t itri = 0; itri < msh.triangles.size(); itri++) {
        /* itri is not only the counter but also the itri number */
        const auto& tri = msh.triangles[itri];
        auto tvtx0 = tri.iv0; /* Index of first vertex */
        auto tvtx1 = tri.iv1; /* Index of second vertex */
        auto tvtx2 = tri.iv2; /* Index of third vertex */
    }

    /**********************************************
     * Data export
     */
    frico::silo db;             // Silo database object
    db.open("topology.silo");   // Open the database
    db.add_mesh("mesh", msh);   // Add the mesh, calling it "mesh"

    std::vector<double> testvar;                // Test variable
    testvar.reserve( msh.triangles.size() );    // For each triangle, compute
    for (auto& tri : msh.triangles) {           // the distance of its barycenter
        auto bar = barycenter(msh, tri);        // from the origin and put
        testvar.push_back( norm(bar) );         // the value in the variable
    }

    // Add the variable to the Silo db, name it "testvar" and associate it
    // to the mesh called "mesh"
    db.add_variable("mesh", "testvar", testvar, frico::var_centering::zonal);

    /**********************************************
     * Graph data structure: here we have a simple data structure
     * to represent graphs. Neither nodes nor edges can be labelled,
     * it just tracks connectivity. It is the graph code that I show
     * to my students, so don't expect too much from it.
     * If we need to extend it, just tell me.
     */
    frico::undirected_graph<size_t> G;

    G.add_edge(1,2);
    G.add_edge(1,3);
    G.add_edge(1,4);
    G.add_edge(1,6);
  
    G.add_edge(2,4);
    G.add_edge(2,5);
    G.add_edge(2,7);

    G.add_edge(3,6);
  
    G.add_edge(4,6);
    G.add_edge(4,7);
    
    G.add_edge(5,7);

    G.add_edge(6,7);

    /* Get a DFS tree, no guarantee about the choice of the root
     * node as the adjacency list is stored in an unordered_map */
    auto dfsT = dfs_tree(G);
    /* Compute the co-tree, subtracting the DFS tree from the initial G */
    auto cotree = G - dfsT;

    /* Save all the three graphs in GraphViz format, to be easily visualized */
    G.save_to_graphviz("graph.dot");
    dfsT.save_to_graphviz("dfs_tree.dot");
    cotree.save_to_graphviz("dfs_cotree.dot");

    /* Compute the fundamental cycles by closing the DFS tree with the
     * co-tree edges taken one by one */
    std::println("Graph cycles:");
    auto cycles = frico::fundamental_cycles(dfsT, cotree);
    for (const auto& cycle : cycles) {
        for (const auto& edge : cycle) {
            std::print(" ({} {}) ", edge.from(), edge.to());
        }
        std::println();
    }

    /* Get a BFS tree, no guarantee about the choice of the root
     * node as the adjacency list is stored in an unordered_map */
    auto bfsT = bfs_tree(G);
    bfsT.save_to_graphviz("bfs_tree.dot");

    return 0;
}