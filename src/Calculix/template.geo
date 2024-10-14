SetFactory("OpenCASCADE");
Point(2) = {1.5, 0, 0, 1.0};
Point(3) = {-1.5, 0, 0, 1.0};
Point(4) = { 3, 2, 0, 1.0};
Point(5) = {-3, 2, 0, 1.0};
Point(6) = {-3, 4, 0, 1.0};
Point(7) = { 0.0, 4, 0, 1.0};
Point(8) = { 3, 4, 0, 1.0};

Point(9) = {-1.5, 2, 0, 1.0};
Point(10) = {1.5, 2, 0, 1.0};
Point(11) = {-0.5, 0, 0, 1.0};
Point(12) = {0.5, 0, 0, 1.0};
Point(13) = {-1, 3.8, 0, 1.0};
Point(14) = {1, 3.8, 0, 1.0};
Point(15) = {0, 3.8, 0, 1.0};

BSpline(1) = {2, 4, 8, 7};
BSpline(2) = {3, 5, 6, 7};
BSpline(5) = {15, 13, 9, 11};
BSpline(6) = {15, 14, 10, 12};

Line(7) = {3, 11};
Line(8) = {11, 12};
Line(9) = {12, 2};

Transfinite Curve {1,2,5,6} = 4*NX Using Progression 1;
Transfinite Curve {7,8,9} = NX Using Progression 1;

Curve Loop(1) = {-2, 1, 9, 6, -5, 7};
Curve Loop(2) = {5, 8, -6};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
