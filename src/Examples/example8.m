%Use GAUSS to compute mass, stiffness matrices and external forces of a tetrahedral mesh
%with a linearly elastic material model, under the influence of gravity
d = [-1 1];
[x,y,z] = meshgrid(d,d,d); % a cube
x = [x(:);0];
y = [y(:);0];
z = [z(:);0];

DT = delaunayTriangulation(x,y,z);

fem = WorldFEM('stvk_linear_tetrahedra', DT.Points, DT.ConnectivityList);

M = mass(fem);
K = stiffness(fem);
f = force(fem);
spy(K);
clear fem
