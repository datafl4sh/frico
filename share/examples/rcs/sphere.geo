//+
SetFactory("OpenCASCADE");
Sphere(1) = {0, -0, 0, 1, -Pi/2, Pi/2, 2*Pi};
MeshSize{ PointsOf{ Surface{:}; } } = 0.1;
