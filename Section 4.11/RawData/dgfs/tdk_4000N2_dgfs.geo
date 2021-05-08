nx1 = 16;
nx2 = 3;
ax1 = 1.0;
ax2 = 1.1;

ny1 = 3;
ny2 = 31;
ny3 = 12;
ay1 = 1.0;
ay2 = 1.0;
ay3 = 1.0;

H0 = 2e-6;

// non-dimensionalized
L = 17e-6 / H0;  // Length of the domain
H = 43e-6 / H0;  // Height of the domain
hG = 2e-6 / H0;  // gap between the wall and the substrate

Lh = 15e-6 / H0;  // length of the heater
Hh = 30e-6 / H0;  // width of the heater

lc = 1; //100;
lcw = 1; //5;
lcb = 1; //12;

Point(1) = {0, 0, 0, lcw};
Point(2) = {Lh, 0, 0, lcw};
Point(3) = {L, 0, 0, lcw};
Point(4) = {L, hG, 0, lcw};
Point(5) = {L, hG+Hh, 0, lcw};
Point(6) = {L, H, 0, lcw};
Point(7) = {Lh, H, 0, lcw};
Point(8) = {0, H, 0, lcw};
Point(9) = {0, hG+Hh, 0, lcw};
Point(10) = {Lh, hG+Hh, 0, lcw};
Point(11) = {Lh, hG, 0, lcw};
Point(12) = {0, hG, 0, lcw};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line(13) = {2, 11};
Line(14) = {11, 4};
Line(15) = {10, 5};
Line(16) = {10, 7};

Line Loop(1) = {1, 13, 11, 12};
Line Loop(2) = {2, 3, -14, -13};
Line Loop(3) = {14, 4, -15, 10};
Line Loop(4) = {15, 5, 6, -16};
Line Loop(5) = {9, 16, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Transfinite Line {1,-11,9,-7} = nx1 Using Progression ax1;
Transfinite Line {2,14,15,-6} = nx2 Using Progression ax2;

Transfinite Line {-12,13,3} = ny1 Using Progression ay1;
Transfinite Line {4,-10} = ny2 Using Progression ay2;
Transfinite Line {5,16,-8} = ny3 Using Progression ay3;

Transfinite Surface {1, 2, 3, 4, 5};
Recombine Surface {1, 2, 3, 4, 5};


Physical Surface("Fluid") = {1, 2, 3, 4, 5};
Physical Line("substrate") = {-11,10,-9};   
Physical Line("symmRight") = {3,4,5};
Physical Line("symmLeftTop") = {-8};    
Physical Line("symmLeftBottom") = {-12};
Physical Line("top") = {6,7};
Physical Line("bottom") = {1,2};

