#include <iostream>

#include "silo.h"

#include "gmsh.h"

#include "geom_mesh.h"
#include "input_gmsh.h"

void
export_mesh(const mommy::mesh& msh)
{
    DBfile *db;
    db = DBCreate("test.silo", DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);

    std::vector<double> x_coords, y_coords, z_coords;
        x_coords.reserve(msh.vertices.size());
        y_coords.reserve(msh.vertices.size());

        for (auto& vtx : msh.vertices) {
            x_coords.push_back( vtx.x() );
            y_coords.push_back( vtx.y() );
            z_coords.push_back( vtx.z() );
        }

        double *coords[] = {x_coords.data(), y_coords.data(), z_coords.data() };

        std::vector<int> nodelist;
        nodelist.reserve( 3*msh.triangles.size() );

        for (auto& t : msh.triangles)
        {
            nodelist.push_back( t.iv0 + 1 );
            nodelist.push_back( t.iv1 + 1 );
            nodelist.push_back( t.iv2 + 1 );
        }

        int lnodelist = nodelist.size();

        int shapetype[] = { DB_ZONETYPE_TRIANGLE };
        int shapesize[] = {3};
        int shapecounts[] = { static_cast<int>(msh.triangles.size()) };
        int nshapetypes = 1;
        int nnodes = msh.vertices.size();
        int nzones = msh.triangles.size();
        int ndims = 3;

        if ( DBPutZonelist2(db, "zonelist", nzones, ndims,
            nodelist.data(), lnodelist, 1, 0, 0, shapetype, shapesize,
            shapecounts, nshapetypes, NULL) < 0 ) {
            std::cout << "DBPutZoneList2() failed" << std::endl;
        }

        if ( DBPutUcdmesh(db, "mesh", ndims, NULL, coords, nnodes,
             nzones, "zonelist", NULL, DB_DOUBLE, NULL) < 0 ) {
            std::cout << "DBPutUcdmesh() failed" << std::endl;
        }

    DBClose(db);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Please specify GMSH .geo file" << std::endl;
        return 1;
    }

    try {
        gmsh::initialize(argc, argv);
        gmsh::open(argv[1]);
    }

    catch (const std::runtime_error& e) {
        std::cerr << "GMSH exception: " << e.what() << std::endl;
        return 1;
    }

    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::setOrder(1);

    mommy::mesh msh;
    mommy::load_mesh_from_gmsh(msh);

    gmsh::finalize();

    std::cout << msh.vertices.size() << " " << msh.triangles.size() << std::endl;

    export_mesh(msh);

    return 0;
}