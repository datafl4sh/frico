SetFactory("OpenCASCADE");
Torus(1) = {0, 0, 0, 0.5, 0.2, 2*Pi};
Torus(2) = {0.95, 0, 0, 0.5, 0.2, 2*Pi};
BooleanUnion{ Volume{1}; Delete;  }{ Volume{2}; Delete; }
MeshSize{ PointsOf{ Volume{:}; } } = 0.05;
Mesh 2;

