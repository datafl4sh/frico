SetFactory("OpenCASCADE");
Torus(1) = {0, -0, 0, 0.5, 0.1, 2*Pi};
MeshSize{ PointsOf{ Volume{1}; } } = 0.1;
Mesh 2;
