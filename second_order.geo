//flow over backward step geometry

Point(1)={0,3,0,1};
Point(2)={10,3,0,1};
Point(3)={10,0,0,1};
Point(4)={17,0,0,1};
Point(5)={17,8,0,1};
Point(6)={0,8,0,1};

//extra points I want to put in close to the walls
//in an attempt at getting some mesh refinement
//along the boundary layers

Point(7)={0,3.1,0,1};
Point(8)={10.1,3.1,0,1};
Point(9)={10.1,0.1,0,1};
Point(10)={17,0.1,0,1};
Point(11)={17,7.9,0,1};
Point(12)={0,7.9,0,1};

//Lines connecting the principal points of the domain
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,1};

//Now form the line loop for the exterior of the domain
Line Loop(7)={1,2,3,4,5,6};

//Create a plane surface that can be meshed
Plane Surface(1)={7};

//Name the surfaces so I can distinguish boundaries and regions
Physical Line(1)={1}; //bottom
Physical Line(2)={2}; //step face
Physical Line(3)={3}; //bottom past step
Physical Line(4)={4}; //outlet
Physical Line(5)={5}; //top
Physical Line(6)={6}; //inlet

//give a name to the whole domain
Physical Surface(10)={1}; //the interior 


Field[1] = MathEval;
Field[1].F = "0.01 + x";
Background Field = 1;
Field[2] = Threshold;
Delete Field [2];
Delete Field [1];
Physical Surface(11) = {1};
Field[1] = MathEval;
Field[1].F = "0.1";
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].Sigmoid = 0;
Background Field = 2;
Field[2].DistMax = 2;
Physical Line(12) = {1, 2};
Background Field = 1;
Background Field = 2;
Delete Field [1];
Field[2].DistMax = 20;
Field[2].DistMin = 10;
Delete Field [2];
Field[1] = Attractor;
Field[1].EdgesList = {1, 2};
Background Field = 1;
Field[2] = Threshold;
Background Field = 2;
Delete Field [2];
Field[2] = MathEval;
Field[2].F = "0.1";
Field[3] = Threshold;
Field[3].DistMax = 5;
Background Field = 3;
Field[3].IField = 1;
