SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
MeshSize{ PointsOf{ Volume{1}; } } = 0.05;
Mesh 2;
