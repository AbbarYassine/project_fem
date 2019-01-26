Mesh.MshFileVersion = 2.2;

R = 1;
N = 10; 
k = 20; //2*Pi; 
h = 2*Pi/(k*N);
xc = 0;
yc = 0;
xc2 = 0;
yc2 = 0;

Point(1) = {xc, yc, 0, h};
Point(2) = {xc + R, yc, 0, h};
Point(3) = {xc, yc + R, 0, h};
Point(4) = {xc - R, yc, 0, h};
Point(5) = {xc, yc - R, 0, h};

Point(6) = {xc2, yc2, 0, h};
Point(7) = {xc2 + 0.5, yc2, 0, h};
Point(8) = {xc2, yc2 + 0.5, 0, h};
Point(9) = {xc2 - 0.5, yc2, 0, h};
Point(10) = {xc2, yc2 - 0.5, 0, h};

Circle(1) = {2,1,3}; // cercle interieur
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {7,6,8}; // cercle interieur
Circle(6) = {8,6,9};
Circle(7) = {9,6,10};
Circle(8) = {10,6,7};

Line Loop(1) = {1:4};       // bord du cercle interieur
Line Loop(2) = {5:8};       // bord du cercle exterieur
Plane Surface(1) = {1,2};   // surface de la couronne
Physical Surface(10) = {1};
Physical Line(11) = {1:4};
Physical Line(12) = {5:8};  