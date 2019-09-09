options ls=80;

data test;
input time x1 $ x2 $;
cards4;
 1 a  A
 2 b  B
10 c  C
50 a  D
 5 b  A
 4 c  B
 8 a  C
40 b  D
60 c  A
20 a  B
21 b  C
22 c  D
 3 a  A
 5 b  B
12 c  C
52 a  D
 7 b  A
 8 c  B
16 a  C
88 a  D
58 a  A
28 a  B
20 a  C
 5 a  B
;;;;

data test2; set test;
   status =1;  * add a dummy status variable;
   if (x2= 'A') then x2a="xA"; else x2a=x2; * A is the reference;
   if (x2= 'B') then x2b="xB"; else x2b=x2; * make B the reference;
   if (x2= 'C') then x2c="xC"; else x2c=x2; * make C the reference;
  
proc phreg data = test2;
   class x1 x2a;
   model time * status(0) = x1 x2a x1*x2a/ type1;

proc phreg data = test2;
   class x1 x2b;
   model time * status(0) = x1 x2b x1*x2b/ type1;

proc phreg data = test2;
   class x1 x2c;
   model time * status(0) = x1 x2c x1*x2c / type1;

proc phreg data = test2;
   class x1 x2d;
   model time * status(0) = x1 x2c x1*x2c / type1;

proc phreg data= test2;
    class x1 x2/ param=effect;
    model time * status(0) = x1 x2 x1*x2 / type1;
