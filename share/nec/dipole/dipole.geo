//+
SetFactory("OpenCASCADE");
Rectangle(1) = {     0,   0,     0,  0.25, 0.01};
Rectangle(2) = { -0.25,   0,     0,  0.25, 0.01};

Coherence;
MeshSize{ PointsOf{ Surface{:}; } } = 0.005;
