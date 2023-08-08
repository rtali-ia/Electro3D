// Gmsh project created on Fri Apr  7 17:13:30 2023
SetFactory("OpenCASCADE");

h = DefineNumber[ 4, Name "h" ];
mxy = DefineNumber[ 20, Name "mxy" ];
mr = DefineNumber[ 30, Name "mr" ];
mz = DefineNumber[ 150, Name "mz" ];

//h = 4; // height of the cylinder
//mxy = 10; // grid in x and y direction
//mr = 15; // grid in radial direction
//mz = 100; // grid in z direction


l1 = Cos(45*Pi/180);
l2 = Cos(22.5*Pi/180);
l3 = Cos(67.5*Pi/180);



//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {l1, l1, 0, 1.0};
//
Point(3) = {l1, -l1, 0, 1.0};
//
Point(4) = {-l1, -l1, 0, 1.0};
//
Point(5) = {-l1, l1, 0, 1.0};
//+
Point(6) = {0.5, 0.5, 0, 1.0};
//
Point(7) = {0.5, -0.5, 0, 1.0};
//
Point(8) = {-0.5, 0.5, 0, 1.0};
//
Point(9) = {-0.5, -0.5, 0, 1.0};
//
Point(10) = {1.0, 0.0, 0, 1.0};
//
Point(11) = {0.0, 1.0, 0, 1.0};
//
Point(12) = {-1.0, 0.0, 0, 1.0};
//
Point(13) = {0.0, -1.0, 0, 1.0};




//
Point(14) = {l2, l3, 0, 1.0};
//
Point(15) = {l3, l2, 0, 1.0};
//
Point(16) = {-l3, l2, 0, 1.0};
//
Point(17) = {-l2, l3, 0, 1.0};
//
Point(18) = {-l2, -l3, 0, 1.0};
//
Point(19) = {-l3, -l2, 0, 1.0};
//
Point(20) = {l3, -l2, 0, 1.0};
//
Point(21) = {l2, -l3, 0, 1.0};


Point(22) = {0.5, 0.25, 0, 1.0};
//
Point(23) = {0.5, -0.25, 0, 1.0};
//
Point(24) = {0.5, 0, 0, 1.0};
//
Point(25) = {-0.5, 0.25, 0, 1.0};
//
Point(26) = {-0.5, -0.25, 0, 1.0};
//
Point(27) = {-0.5, 0, 0, 1.0};
//
Point(28) = {0.25, -0.5, 0, 1.0};
//
Point(29) = {-0.25, -0.5 , 0, 1.0};
//
Point(30) = {0, -0.5, 0, 1.0};
//
Point(31) = {0.25, 0.5, 0, 1.0};
//
Point(32) = {-0.25, 0.5 , 0, 1.0};
//
Point(33) = {0, 0.5, 0, 1.0};



Point(34) = {0.25, 0.25, 0, 1.0};
//
Point(35) = {-0.25, -0.25, 0, 1.0};
//
Point(36) = {-0.25, 0.25 , 0, 1.0};
//
Point(37) = {0.25, -0.25, 0, 1.0};
//
Point(38) = {0.25, 0.0, 0, 1.0};
//
Point(39) = {-0.25, 0.0, 0, 1.0};
//
Point(40) = {0.0, 0.25 , 0, 1.0};
//
Point(41) = {0.0, -0.25, 0, 1.0};


//+
Circle(1) = {17, 1, 5};
//+
Circle(2) = {5, 1, 16};
//+
Circle(3) = {16, 1, 11};
//+
Circle(4) = {11, 1, 15};
//+
Circle(5) = {15, 1, 2};
//+
Circle(6) = {2, 1, 14};
//+
Circle(7) = {14, 1, 10};
//+
Circle(8) = {10, 1, 21};
//+
Circle(9) = {21, 1, 3};
//+
Circle(10) = {3, 1, 20};
//+
Circle(11) = {20, 1, 13};
//+
Circle(12) = {13, 1, 19};
//+
Circle(13) = {19, 1, 4};
//+
Circle(14) = {4, 1, 18};
//+
Circle(15) = {18, 1, 12};
//+
Circle(16) = {12, 1, 17};
//+
Line(17) = {9, 26};
//+
Line(18) = {26, 27};
//+
Line(19) = {27, 25};
//+
Line(20) = {25, 8};
//+
Line(21) = {8, 32};
//+
Line(22) = {32, 36};
//+
Line(23) = {36, 25};
//+
Line(24) = {36, 40};
//+
Line(25) = {40, 34};
//+
Line(26) = {34, 22};
//+
Line(27) = {22, 6};
//+
Line(28) = {6, 31};
//+
Line(29) = {31, 33};
//+
Line(30) = {33, 32};
//+
Line(31) = {40, 33};
//+
Line(32) = {40, 1};
//+
Line(33) = {1, 41};
//+
Line(34) = {41, 30};
//+
Line(35) = {36, 39};
//+
Line(36) = {39, 35};
//+
Line(37) = {35, 29};
//+
Line(38) = {31, 34};
//+
Line(39) = {34, 38};
//+
Line(40) = {38, 37};
//+
Line(41) = {37, 28};
//+
Line(42) = {28, 7};
//+
Line(43) = {7, 23};
//+
Line(44) = {23, 24};
//+
Line(45) = {24, 22};
//+
Line(46) = {24, 38};
//+
Line(47) = {38, 1};
//+
Line(48) = {1, 39};
//+
Line(49) = {39, 27};
//+
Line(50) = {26, 35};
//+
Line(51) = {35, 41};
//+
Line(52) = {41, 37};
//+
Line(53) = {37, 23};
//+
Line(54) = {28, 30};
//+
Line(55) = {30, 29};
//+
Line(56) = {29, 9};
//+
Line(57) = {8, 5};
//+
Line(58) = {32, 16};
//+
Line(59) = {33, 11};
//+
Line(60) = {31, 15};
//+
Line(61) = {6, 2};
//+
Line(62) = {14, 22};
//+
Line(63) = {24, 10};
//+
Line(64) = {23, 21};
//+
Line(65) = {7, 3};
//+
Line(66) = {20, 28};
//+
Line(67) = {30, 13};
//+
Line(68) = {29, 19};
//+
Line(69) = {4, 9};
//+
Line(70) = {18, 26};
//+
Line(71) = {27, 12};
//+
Line(72) = {25, 17};
//+
Curve Loop(1) = {3, -59, 30, 58};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, -60, 29, 59};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, -61, 28, 60};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6, 62, 27, 61};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {7, -63, 45, -62};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {8, -64, 44, 63};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {9, -65, 43, 64};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {10, 66, 42, 65};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {11, -67, -54, -66};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {55, 68, -12, -67};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {56, -69, -13, -68};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {70, -17, -69, 14};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {71, -15, 70, 18};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {72, -16, -71, 19};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {57, -1, -72, 20};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {58, -2, -57, 21};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {22, 23, 20, 21};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {35, 49, 19, -23};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {36, -50, 18, -49};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {37, 56, 17, 50};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {37, -55, -34, -51};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {52, 41, 54, -34};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {53, -43, -42, -41};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {46, 40, 53, 44};
//+
Plane Surface(24) = {24};
//+
Curve Loop(25) = {47, 33, 52, -40};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {48, 36, 51, -33};
//+
Plane Surface(26) = {26};
//+
Curve Loop(27) = {24, 32, 48, -35};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {25, 39, 47, -32};
//+
Plane Surface(28) = {28};
//+
Curve Loop(29) = {26, -45, 46, -39};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {28, 38, 26, 27};
//+
Plane Surface(30) = {30};
//+
Curve Loop(31) = {29, -31, 25, -38};
//+
Plane Surface(31) = {31};
//+
Curve Loop(32) = {30, 22, 24, 31};
//+
Plane Surface(32) = {32};
//+
Transfinite Curve {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 21, 20, 23, 22, 35, 19, 49, 50, 36, 18, 17, 56, 37, 55, 34, 51, 52, 41, 54, 42, 43, 53, 40, 33, 48, 32, 47, 39, 46, 45, 44, 27, 26, 38, 25, 31, 24, 30, 29, 28} = mxy Using Progression 1;
//+
Transfinite Curve {59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 57, 58} = mr Using Progression 1;
//+
Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
//+
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};


newEntities[]=
Extrude {0, 0, h}
{
	Surface{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
	Layers{mz};
	Recombine;
};
//+
Physical Surface("D1", 1) = {50};
//+
Physical Surface("D2", 2) = {42};
//+
Physical Surface("D3", 3) = {33};
//+
Physical Surface("D4", 4) = {91};
//+
Physical Surface("D5", 5) = {83};
//+
Physical Surface("D6", 6) = {76};
//+
Physical Surface("D7", 7) = {66};
//+
Physical Surface("D8", 8) = {58};
//+
Physical Surface("N", 9) = {54, 62, 72, 80, 87, 94, 38, 46};
//+
Physical Surface("top", 10) = {45, 41, 37, 96, 93, 89, 85, 81, 99, 136, 135, 133, 49, 53, 57, 61, 65, 115, 118, 131, 129, 126, 102, 69, 73, 77, 113, 121, 110, 123, 105, 107};
//+
Physical Surface("bottom", 11) = {5, 4, 3, 2, 1, 16, 15, 32, 31, 30, 29, 24, 23, 22, 28, 17, 25, 26, 14, 13, 12, 18, 27, 19, 20, 21, 6, 7, 8, 9, 10, 11};
//+
Physical Volume("interior", 12) = {3, 4, 5, 30, 2, 1, 32, 31, 29, 25, 28, 27, 18, 17, 16, 24, 6, 7, 8, 23, 22, 21, 26, 19, 20, 10, 9, 11, 12, 13, 14, 15};
