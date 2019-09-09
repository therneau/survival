libname  save   '/bigdata/users/therneau/sasdata';
filename whs1   '/bigdata/users/therneau/sasdata/nhs1.csv';
filename whs2   '/bigdata/users/therneau/sasdata/nhs2.csv';


proc import datafile=whs1
    out = save.nhs
    dbms= csv  replace;

proc import datafile=whs2
    out= save.nhs2 dbms=csv replace;
