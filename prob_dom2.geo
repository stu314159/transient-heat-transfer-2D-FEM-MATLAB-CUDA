//define a problem domain for the gmsh program:
//basic geometry will be a rectangle with a circlular hole in the middle
//to be used with my transient heat transfer fem code

//these are the four points that define the rectangular domain
Point(1)={0,0,0,.1};
Point(2)={5,0,0,.1};
Point(3)={5,3,0,.1};
Point(4)={0,3,0,.1};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

//the circle will be at {2.5,1.5,0} the center of the rectangle

Point(5)={2.5,1.5,0,.05}; //center point
Point(6)={3,1.5,0,.05}; //right point
Point(7)={2.5,2,0,.05}; //top point
Point(8)={2,1.5,0,.05}; //left point
Point(9)={2.5,1,0,.05}; //bottom point

//now form the Circle segments:
Circle(5)={6,5,7};
Circle(6)={7,5,8};
Circle(7)={8,5,9};
Circle(8)={9,5,6};

//Now form the line loops for the exterior of the rectangle
Line Loop(9)={1,2,3,4};

//...and the exterior of the circular hole
Line Loop(10)={5,6,7,8};

//Create a plane surface that can be meshed:
Plane Surface(11)={9,10};

//Name my sufurfaces so that I can distinguish boundaries and regions
Physical Line(12) = {1}; //bottom
Physical Line(13) = {2}; //right
Physical Line(14) = {3}; //top
Physical Line(15) = {4}; //left

//I wonder if I can just assign a single number to a line loop
Physical Line(16)={5};
Physical Line(17)={6};
Physical Line(18)={7};
Physical Line(19)={8};

//ok, now I need to give a name to the whole prroblem domain.
Physical Surface(20)={11};



