# This generates a mesh to replicate the results of the
# first numerical experiment described in the section 13
# of the MOM3D manual

SetFactory("OpenCASCADE");

c = 299792458;
f = 3.0e9;
lambda = c/f;
w = 1.4*lambda;
h = 1.4*lambda;

Rectangle(1) = {-w/2, -h/2, 0, w, h, 0};

MeshSize{ PointsOf{ Surface{:}; } } = 0.005;
