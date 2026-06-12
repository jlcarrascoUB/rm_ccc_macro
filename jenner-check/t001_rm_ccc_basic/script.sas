/* ---------------------------------------------------------------
   Caller for %uccc — U-Statistics Concordance Correlation
   Coefficient for repeated longitudinal measures.
   (rm_ccc_macro_v2.sas, Carrasco et al. 2013)

   Demonstrates the %uccc macro from rm_ccc_macro_v2.sas using
   mock blood-pressure readings: 16 subjects, 2 raters (methods),
   3 time points, identity weight matrix.
   --------------------------------------------------------------- */

options ls=80 ps=65;

/* --- %ic_ccc helper (from rm_ccc_macro_v2.sas) --- */
%macro ic_ccc(CCC,dev,S);
  SECCC = sqrt(&dev*&S*&dev`);
  Z=0.5*log((1+&CCC)/(1-&CCC));
  SEZ=sqrt((SECCC**2)/(((1+&CCC)**2)*((1-&CCC)**2)));
  T=1.95996; /* probit(0.975) */
  LLZ=Z-T*SEZ;
  ULZ=Z+T*SEZ;
  LL95 = (exp(2*LLZ)-1)/(exp(2*LLZ)+1);
  UL95 = (exp(2*ULZ)-1)/(exp(2*ULZ)+1);
  result= &CCC || LL95 || UL95 || SECCC || Z || SEZ;
%mend;

/* --- %uccc macro (from rm_ccc_macro_v2.sas, identity-weight path) --- */
%MACRO UCCC(DATA, VAR1, VAR2, VAR3, weight=, type=, output=uout, delta=);
proc sort data=&data out=_temp; by &var3; run;
proc transpose data=_temp out=_temp2 prefix=temp;
  var &var1;
  by &var3;
run;
proc transpose data=_temp2 out=_x prefix=x; var temp:; run;

proc transpose data=_temp out=_temp2 prefix=temp;
  var &var2;
  by &var3;
run;
proc transpose data=_temp2 out=_y prefix=y; var temp:; run;

proc iml;
  use _x; read all var _num_ into X;
  use _y; read all var _num_ into Y;
  xy = x||y;
  N  = NROW(XY);
  p  = ncol(xy)/2;
  labels={"CCC" "Lower 95% CL" "Upper 95% CL" "SE CCC" "Z" "SE Z"};

  use &weight; read all var _num_ into &type;
  print 'Repeated Measures Concordance Correlation Coefficient',
        'Derived from U-Statistics';

  PSI1_I = J(N,1,0);
  PSI2_I = J(N,1,0);
  DO I = 1 TO N;
    PSI1 = 0; PSI2 = 0;
    DO J = 1 TO N;
      B = 0.5*( (abs(X[J,1:p]-Y[J,1:p])##&delta)*D1*
                (abs(X[J,1:p]-Y[J,1:p])##&delta)` );
      A = 0.5*( (abs(X[I,1:p]-Y[J,1:p])##&delta)*D1*
                (abs(X[I,1:p]-Y[J,1:p])##&delta)` ) +
          0.5*( (abs(X[J,1:p]-Y[I,1:p])##&delta)*D1*
                (abs(X[J,1:p]-Y[I,1:p])##&delta)` );
      PSI1 = PSI1 + B;
      PSI2 = PSI2 + A;
    END;
    BI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D1*
         (abs(X[I,1:p]-Y[I,1:p])##&delta)`;
    AI = BI;
    PSI1_I[I] = PSI1 + ((N-2)/2)#BI;
    PSI2_I[I] = PSI2 - AI;
  END;

  U       = SUM(PSI1_I)/(N*(N-1));
  V       = SUM(PSI2_I)/(N*(N-1));
  RHOCRM1 = ((N-1)*(V-U))/(U+(N-1)*V);

  PSI1I = PSI1_I/(N-1);
  PSI2I = PSI2_I/(N-1);
  PSI_I = (PSI1I||PSI2I)`;
  Uvec  = U//V;
  U_M   = REPEAT(Uvec,1,N);
  DELTA = PSI_I - U_M;
  c     = {1 2, 2 4};
  V_U   = (c/(N)##2)#(DELTA*DELTA`);
  partialU = (-n#(n-1)#V)/(U+(n-1)#V)##2;
  partialV = ( n#(n-1)#U)/(U+(n-1)#V)##2;
  d     = partialU||partialV;
  V_RHO = d*V_U*d`;
  SE_RHO= SQRT(V_rho);
  T     = 1.95996;  /* probit(0.975) */
  Z     = LOG((1+rhocrm1)/(1-rhocrm1))/2;
  SE_Z  = SQRT(V_rho)/(1-rhocrm1##2);
  Z_LCL = Z-(SE_Z#T);
  Z_UCL = Z+(SE_Z#T);
  LCL1  = (EXP(2#Z_LCL)-1)/(EXP(2#Z_LCL)+1);
  UCL1  = (EXP(2#Z_UCL)-1)/(EXP(2#Z_UCL)+1);

  output = rhocrm1 || lcl1 || ucl1 || se_rho || Z || se_z;
  print output[label='' colname=labels rowname=' ' format=15.5];
  print &type[label='Weight Matrix'];

  create &output from output[colname=labels];
  append from output;
quit;
%MEND UCCC;

/* --- Identity weight matrix (3 x 3 for 3 time points) --- */
data _d;
  col1=1; col2=0; col3=0; output;
  col1=0; col2=1; col3=0; output;
  col1=0; col2=0; col3=1; output;
run;

/* --- Mock data: 16 subjects, 2 methods, 3 time points --- */
/* One row per (subject, method, time); y = measurement      */
data meas;
  input subj method time y;
  datalines;
1 1 1 120.1
1 1 2 121.4
1 1 3 119.8
1 2 1 121.3
1 2 2 122.5
1 2 3 121.0
2 1 1 130.5
2 1 2 131.2
2 1 3 129.9
2 2 1 131.8
2 2 2 132.4
2 2 3 131.1
3 1 1 115.0
3 1 2 116.3
3 1 3 114.7
3 2 1 116.2
3 2 2 117.4
3 2 3 115.9
4 1 1 125.3
4 1 2 126.1
4 1 3 124.8
4 2 1 126.5
4 2 2 127.3
4 2 3 126.0
5 1 1 118.7
5 1 2 119.5
5 1 3 118.2
5 2 1 119.9
5 2 2 120.7
5 2 3 119.4
6 1 1 135.2
6 1 2 136.0
6 1 3 134.8
6 2 1 136.4
6 2 2 137.2
6 2 3 136.0
7 1 1 122.0
7 1 2 122.8
7 1 3 121.5
7 2 1 123.2
7 2 2 124.0
7 2 3 122.7
8 1 1 128.4
8 1 2 129.2
8 1 3 127.9
8 2 1 129.6
8 2 2 130.4
8 2 3 129.1
9 1 1 117.6
9 1 2 118.4
9 1 3 117.1
9 2 1 118.8
9 2 2 119.6
9 2 3 118.3
10 1 1 132.9
10 1 2 133.7
10 1 3 132.4
10 2 1 134.1
10 2 2 134.9
10 2 3 133.6
11 1 1 119.3
11 1 2 120.1
11 1 3 118.8
11 2 1 120.5
11 2 2 121.3
11 2 3 120.0
12 1 1 127.7
12 1 2 128.5
12 1 3 127.2
12 2 1 128.9
12 2 2 129.7
12 2 3 128.4
13 1 1 113.4
13 1 2 114.2
13 1 3 112.9
13 2 1 114.6
13 2 2 115.4
13 2 3 114.1
14 1 1 124.8
14 1 2 125.6
14 1 3 124.3
14 2 1 126.0
14 2 2 126.8
14 2 3 125.5
15 1 1 131.6
15 1 2 132.4
15 1 3 131.1
15 2 1 132.8
15 2 2 133.6
15 2 3 132.3
16 1 1 116.9
16 1 2 117.7
16 1 3 116.4
16 2 1 118.1
16 2 2 118.9
16 2 3 117.6
;
run;

/* Pivot to wide form: one row per (subject, time), columns met1 met2 */
proc sort data=meas; by subj time method; run;
proc transpose data=meas out=meas_wide prefix=met;
  by subj time;
  var y;
  id method;
run;

/* --- Run %uccc (U-statistics CCC, identity weight matrix) --- */
%uccc(meas_wide, met1, met2, time,
      weight=_d, type=d1, output=ccc_out, delta=1);

/* --- Print results --- */
title "Concordance Correlation Coefficient - Repeated Measures";
title2 "U-Statistics method, identity weight (Carrasco et al. 2013)";
proc print data=ccc_out noobs label; run;
title; title2;
