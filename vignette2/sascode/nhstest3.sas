libname save "/bigdata/users/therneau/sasdata";
options ls=80;

proc phreg data= save.nhs;
    where id < 40566;     * 1/3 the subjects;
    model (age1 age2) * death(0) = x1 x2 xb_1 xb_2 xb_3 xb_4 xb_5 xb_6 xb_7
        xb_8 xb_9 xb_10 xb_11 xb_12 xb_13 xb_14 xb_15 xb_16 xb_17 xb_18 xb_19
        xb_20 xb_21 xb_22 xb_23 xb_24 xb_25 xb_26 xb_27
        xb_28 xb_29 xb_30 xb_31 xb_32 xb_33 xb_34 xb_35 xb_36
        xb_37 xb_38 xb_39 xb_40;
    strata agem;

proc phreg data= save.nhs;
    where id < 81133;   * 2/3;    
    model (age1 age2) * death(0) = x1 x2 xb_1 xb_2 xb_3 xb_4 xb_5 xb_6 xb_7
        xb_8 xb_9 xb_10 xb_11 xb_12 xb_13 xb_14 xb_15 xb_16 xb_17 xb_18 xb_19
        xb_20 xb_21 xb_22 xb_23 xb_24 xb_25 xb_26 xb_27
        xb_28 xb_29 xb_30 xb_31 xb_32 xb_33 xb_34 xb_35 xb_36
        xb_37 xb_38 xb_39 xb_40;
    strata agem;
    
proc phreg data= save.nhs;
    model (age1 age2) * death(0) = x1 x2 xb_1 xb_2 xb_3 xb_4 xb_5 xb_6 xb_7
        xb_8 xb_9 xb_10 xb_11 xb_12 xb_13 xb_14 xb_15 xb_16 xb_17 xb_18 xb_19
        xb_20 xb_21 xb_22 xb_23 xb_24 xb_25 xb_26 xb_27
        xb_28 xb_29 xb_30 xb_31 xb_32 xb_33 xb_34 xb_35 xb_36
        xb_37 xb_38 xb_39 xb_40;
    strata agem;

proc phreg data=save.nhs;
    model (age1 age2)*death(0) = x1 x2 xb_1 xb_2 xb_3 xb_4 xb_5 xb_6 xb_7
        xb_8 xb_9 xb_10 xb_11 xb_12 xb_13 xb_14 xb_15 xb_16 xb_17 xb_18 xb_19
        xb_20 xb_21 xb_22 xb_23 xb_24 xb_25 xb_26 xb_27
        xb_28 xb_29 xb_30 xb_31 xb_32 xb_33 xb_34 xb_35 xb_36
        xb_37 xb_38 xb_39 xb_40;
    strata month1;

proc phreg data=save.nhs;
    model age2*death(0) = x1 x2 xb_1 xb_2 xb_3 xb_4 xb_5 xb_6 xb_7
        xb_8 xb_9 xb_10 xb_11 xb_12 xb_13 xb_14 xb_15 xb_16 xb_17 xb_18 xb_19
        xb_20 xb_21 xb_22 xb_23 xb_24 xb_25 xb_26 xb_27
        xb_28 xb_29 xb_30 xb_31 xb_32 xb_33 xb_34 xb_35 xb_36
        xb_37 xb_38 xb_39 xb_40;
    strata month1;

proc phreg data=save.nhs;
    model month2*death(0) = x1 x2 xb_1 xb_2 xb_3 xb_4 xb_5 xb_6 xb_7
        xb_8 xb_9 xb_10 xb_11 xb_12 xb_13 xb_14 xb_15 xb_16 xb_17 xb_18 xb_19
        xb_20 xb_21 xb_22 xb_23 xb_24 xb_25 xb_26 xb_27
        xb_28 xb_29 xb_30 xb_31 xb_32 xb_33 xb_34 xb_35 xb_36
        xb_37 xb_38 xb_39 xb_40;
    strata month1;

