//+
SetFactory("OpenCASCADE");

width = 0.002*3.141592;

Rectangle(1) = {     0,   0,     0,  0.25, width};
Rectangle(2) = { -0.25,   0,     0,  0.25, width};

Coherence;
MeshSize{ PointsOf{ Surface{:}; } } = 0.005;
