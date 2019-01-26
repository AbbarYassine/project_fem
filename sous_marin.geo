
DefineConstant[
 h = {0.05, Min 0.01, Max 10, Step 0.01, Name "h"} 
];

//Garder finesse du maillage

Mesh.CharacteristicLengthMax=h;
Mesh.CharacteristicLengthMin=h;

//Ellipse
Point(1)= {1,1,0,h}; 
Point(2)= {2,1,0,h}; 
Point(3)= {1,1.8,0,h}; 
Point(4)= {0,1,0,h}; 
Point(5)={1,0.2,0,h};

Ellipse(1) = {2,1,2,3}; 
Ellipse(2) = {3,1,2,4}; 
Ellipse(3) = {4,1,2,5}; 
Ellipse(4) = {5,1,2,2};

Line Loop(1) = {2,3,4,1}; //Pourtour ellipse

//Sous-marin

Point(6)= {0.7,1.1,0,h}; 
Point(7)= {1.1,1.1,0,h}; 

Point(8)= {1.2,1.1,0,h};
Point(9)= {1.4,1.1,0,h}; 

Point(10)={0.7,0.9,0,h}; 
Point(11)={1.4,0.9,0,h};

Point(12)={1.5,1,0,h};
Point(13)={1.4,1,0,h};//centre c103 et c104

Point(14)={0.5,1,0,h}; 

Line(100)={6,7}; 
Line(101)={8,9}; 
Line(102)={10,11}; 

Circle(103) = {12,13,9};
Circle(104) = {12,13,11};

//Tete
Point(15)={0.65,1.2,0,h};
Line(105)={6,15}; 

Point(16) = {0.6,1.15,0,h}; 
Line(106)={15,16}; 

Point(17) = {0.6,1.1,0,h}; 
Line(107)={16,17};

Point(18) = {0.65,0.8,0,h}; 
Line(108)={10,18};

Point(19) = {0.6,0.85,0,h}; 
Line(109)={18,19};

Point(20) = {0.6,0.9,0,h}; 
Line(110)={19,20};

Point(21) = {0.7,1,0,h};

Line(111)={14,20};

Point(23) = {0.6,1,0,h}; 
Circle(112) = {14,23,17};
Circle(113) = {17,23,21};
Circle(114) = {21,23,20};

Line(115)={14,17};
Line(116)={17,21};
Line(117)={21,20};

Circle(118) = {20,23,14};

Point(24) = {1.05,1.2,0,h};
Line(119)={7,24};

Point(25) = {1.15,1.17,0,h};
Line(120)={8,25};

Point(26) = {1.1,1.185,0,h}; 
Line(121) = {24,26}; 
Line(122) = {26,25};

