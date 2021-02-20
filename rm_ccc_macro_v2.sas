/********************************************************************************

Computation of a Concordance Correlation Coefficient for longitudinal repeated
measurements via U-statisitics and/or a variance components linear mixed model. 

*********************************************************************************

MACRO DESCRIPTION:															
--------------------------------------------------------------------------------- 	
This macro produces a concordance correlation coefficient (CCC) for longitudinal 
repeated measurements. The CCC computed from both (1) U-statistic methodology 
(King et al., 2007) and (2) through an appropriate specification of the intra-
class correlation coefficient from a variance components linear mixed model 
(Carrasco et al., 2009). Both are produced, except in some specific situations as 
follows.  For missing or unbalanced data, or for multiple observers/raters, the 
U-statistic method is not appropriate.  For small samples (n<20), both CCC's are 
produced, but the U-statistics method is recommended, and this is noted in the 
log.  (The U-statistics method is also recommended for non-normal or qualitative 
data, though this is left to the user's discretion.)  Both methods are approp-
riate for both an identity or a diagonal weight matrix, D, but because the VC 
approach is not currently developed for a general matrix, only the U-statistics 
approach is valid if the D matrix is of a general form. If repeated measures are 
replicated, but not longitudinal in nature, then this is indicated by the NL 
option as described below to produce the correct result. 

The macro does not allow for covariates at this time.

DATASET REQUIREMENTS:															
--------------------------------------------------------------------------------- 	
The dataset must have one line per measurement, and that line should contain
the measurement variable, observer/method variable, a discrete time variable, and 
subject identification variable.  

MACRO VARIABLES:																
--------------------------------------------------------------------------------- 	
DATA=			Input data set, as described above.

Y=				Measurement variable 

METHOD=			Method, or observer, variable 

TIME=			Discrete time variable, or in the case of non-longitudinal data, 
				a variable differentiating the repeated measures. 

SUBJECT=		Subject identifier 

WEIGHT=			A data set with p rows and p columns representing a pxp general 
				weight matrix, where p is the number of repeated measures per 
				method and subject.  If no weight matrix is specified, and no 
				diagonal weight matrix is specified, then an identity matrix is 
				used. (See also DIAG.)

DIAG=			An alternative way to enter a weight matrix of diagonal form. A 
				list of p weights should be listed within parentheses, separated 
				by commas (where p is the number of repeated measures per method 
				and subject).  If no weight matrix is specified through either 
				the WEIGHT or DIAG options, then an identity matrix is used.  
				(See also WEIGHT.)

				For example, to specify a diagonal weight matrix for a dataset 
				with 3 repeated measures per method and subject, the following 
				macro call could be used: 
				%rm_ccc(data=data.example,y=y,method=met,time=t,subject=subj,diag=(1,2,3)); 

NL=				0 if data is longitudinal, 1 if data is nonlongitudinal.  
				The default is 0.

TYPE=			TYPE=0 for no repeated statement, TYPE=1 for autoregressive 
				(AR(1)), or TYPE=2 for compound symmetry covariance structure.  
				Note that this option has no effect if NL=1.

OUTPUT=			Specifies an output dataset.  If no dataset is specified, no 
				output dataset is created.   

MODEL_DETAIL=	OFF for simple output.  If ON is specified, then model details 
				will be printed for the linear mixed model used in the variance 
				component method.

DELTA=			

********************************************************************************/

%macro rm_ccc(data=,y=,method=,time=,subject=,weight=,diag=,nl=0, type=0,  
			    output=, model_detail=OFF, delta=1);

data _calc_ccc;
  set &data;
  keep &y &method &time &subject;
run;

%put model_detail is: &model_detail;

*Code for identity matrix as the default*;
%let u=d1; %let vc=vc_unw8ed;

proc sort data=_calc_ccc;
  by &subject &time &method;
run;
proc sql noprint;
  select count(distinct &method), nmiss(&y)
    into :nmethods, :nmiss
	  from _calc_ccc;
  select count(distinct &subject) into :n from _calc_ccc;
  create table _balanced as 
    select distinct &subject, count(&y) as ni, count(distinct &time) as p, 
		   count(distinct &method) as k
	  from _calc_ccc
	    group by &subject;
  select max(p), min(p), max(k), min(k), max(ni), min(ni)
    into :maxp, :minp, :maxk, :mink, :maxn, :minn
	  from _balanced;
quit;

/* Check properties of dataset to determine proper methodology */
%if &nmethods>2 |&nmiss>0 | &maxp ne &minp | &maxk ne &mink | 
    &maxn ne &minn %then %do;
	%put;
	%put NOTE: THE U-STATISTICS METHODOLOGY IS NOT YET DEVELOPED FOR MISSING OR UNBALANCED DATA, OR FOR MORE THAN 2 METHODS.;
    %put;
	%let u=0; 
	%if &n <20 %then %do;
	  %put;
	  %put WARNING: VARIANCE COMPONENTS METHODOLOGY MAY NOT BE RECOMMENDED DUE TO SMALL SAMPLE SIZE (N<20).;
	  %put;
	%end;
%end;
%else %if &n < 20 %then %do;
	%put;
	%put NOTE: THE U-STATISTIC METHODOLOGY IS RECOMMENDED DUE TO SMALL SAMPLE SIZE (N<20).;
	%put;
%end;

/* Check for correct format of diag macro variable input */ 
%if &diag ne %then %do;
  %let noparens=%scan(&diag,1,%nrstr(()) );
  %if %quote(&diag)=%quote(&noparens) %then %do;
	%put;
	%put WARNING: THE DIAGONAL MATRIX IS OF THE WRONG FORMAT. AN IDENTITY MATRIX WILL BE USED.;
	%put;
    %let diag= ;
  %end;
%end;

/* Check weight matrix properties */
%if &nl=1 %then %do;
  %let weight=; %let diag=;
%end;
%if &diag ne %then %do; %let weight= ; %end;

proc iml;
  *identity matrix is the default*;
  d=I(&maxp);
  
  /********** Set up macro variables for NONLONGITUDINAL data ****************/
  %if &nl=1 %then %do; 
    %if &u ^= 0 %then %let u=d4; 
	%let vc=ccc;
	d=J(&maxp,&maxp,1);
  %end;
  /********** Set up macro variables for a DIAGONAL weight matrix ************/
  %if &diag ne %then %do; 
    d=block&diag;
    ncol=ncol(d);
	call symputx('ncol',ncol);
    %if &ncol ^= &maxp %then %do;
	  %put;
	  %put ERROR: THE DIAGONAL MATRIX IS OF THE WRONG DIMENSION. THE MACRO WILL TERMINATE;
	  %put;
	  %let u=0; %let vc=0;
    %end;
	%else %do;
	  %if &u ^= 0 %then %let u=d2; %let vc=vc_w8ed;
	%end;
  %end;
  /****** Set up macro variables for a WEIGHT matrix (diagonal or general) ****/
  %if &weight ne %then %do; 
	use &weight; read all into d;
	if d=diag(d) then call symputx('indic',1);
	  else call symputx('indic',0);
	%if &indic %then %do;
	  %if &u ^= 0 %then %let u=d2; %let vc=vc_w8ed;
  	%end;
	%else %do;
	  %put;
	  %put NOTE: VARIANCE COMPONENTS METHODOLOGY IS NOT YET DEVELOPED FOR A GENERAL WEIGHT MATRIX.;
	  %put;
	  %if &u ^= 0 %then %let u=d3; %let vc=0; 
    %end;
  %end;

  ncol=ncol(d);
  call symputx('ncol',ncol);
  %if &ncol ^= &maxp %then %do;
  	%put;
    %put ERROR: THE WEIGHT MATRIX OF THE WRONG DIMENSION. THE MACRO WILL TERMINATE.; 
	%put;
    %let u=0; %let vc=0;
  %end;
  *********************************************************************************;

  create _d from d; append from d;
quit;

%if (&type <0 or &type >2) %then %do;
	%put;
	%put ERROR: TYPE MUST BE 0, 1 OR 2;
	%put;
	%let vc=0;
%end;

%if &u=0 and &vc=0 %then %do;
	%put;
	%put ERROR: NEITHER VARIANCE COMPONENTS NOR U-STATISTICS METHODOLOGY IS VALID. PROGRAM WILL TERMINATE;
	%put;
%end;

*** Calculate CCC via U-statistics methodology ***;
%if &u ^= 0 %then %do;
  proc transpose data=_calc_ccc out=_u prefix=met;
    by &subject &time;
    var &y;
  run;
  %uccc(_u,met1,met2,&time,weight=_d,type=&u,output=_uout,delta=&delta);
%end;

*** Calculate CCC via Variance Components methodology ***;
%if &vc ne 0 %then %do;
  data _vc;
    set _calc_ccc;
    y=&y;
    met=&method;
    time=&time;
    ind=&subject;
  run;

  /***** Unweighted variance components method *****/
  %if &vc=vc_unw8ed %then %do;
    %ccc_lon(_vc,_vcout,rho=&type,D=0,detail=&model_detail); 
  %end;
  /***** Weighted variance components method ******/
  %else  %if &vc=vc_w8ed %then %do;
    %ccc_lon(_vc,_vcout,rho=&type,D=_d,detail=&model_detail);
  %end;
  /***** Nonlongitudinal variance components method *****/
  %else %if &vc=ccc %then %do;
	%ccc(_vc,output=_vcout,m=&maxp,detail=&model_detail);
  %end;
%end;
ods select none;
proc datasets;  
  delete _x _y _vc _u _temp2 _temp _table1 _table2 _table3 _table4 _table5 _table6 _d _balanced; 
quit;
ods select all;
%if &output ne %then %do;
  %if &u ne 0 and &vc ne 0 %then %do; data &output; merge _uout _vcout; run;  %end;
  %else %if &u ne 0 %then %do; data &output; set _uout; run; %end;
  %else %if &vc ne 0 %then %do; data &output; set _vcout; run; %end;
%end;

%if &output ne _vcout and &output ne _uout %then %do; proc datasets; delete _vcout _uout; quit; %end;

%mend rm_ccc;
 


OPTIONS LS=80 PS=65;
************************************************************************************
%UCCC produces the CCC for Identity, Diagonal, and General Covariance Matrices, as 
well as for Non-Longitudinal Repeated Measures data 
************************************************************************************;
%MACRO UCCC(DATA, VAR1, VAR2, VAR3, weight=, type=, output=uout,delta=);
proc sort data=&data out=_temp; by &var3; run;
proc transpose data=_temp out=_temp2 prefix=temp;
  var &var1;
  by &var3;
run;
proc transpose data=_temp2 out=_x prefix=x; var temp:;

proc transpose data=_temp out=_temp2 prefix=temp;
  var &var2;
  by &var3;
run;
proc transpose data=_temp2 out=_y prefix=y; var temp:;

proc iml;
use _x;
read all into X;
use _y;
read all into Y;
xy = x||y;
N = NROW(XY);
p=ncol(xy)/2;
mux = (x[+,]/n)`;
muy = (y[+,]/n)`;

sum_xy=xy[+,];
*dividing by n here instead of (n-1) so that this approach and U-stat approach agree;
/*cov_est=(xy`*xy-sum_xy`*sum_xy/n)/(n); 
xxhat = cov_est[1:p,1:p];
yyhat = cov_est[p+1:2*p,p+1:2*p];
xyhat = cov_est[1:p,p+1:2*p];
yxhat = cov_est[p+1:2*p,1:p];
*/
use &weight;
read all into &type;
/*num = trace(&type*xyhat + &type*yxhat);
denom = trace(&type*xxhat + &type*yyhat)+ (mux-muy)`*&type*(mux-muy);
rho = num/denom;
*/
print 'Repeated Measures Concordance Correlation Coefficient Derived from U-Statistics';*,,
        rho[label='' colname=('CCC')];
labels={"CCC" "Lower 95% CL" "Upper 95% CL" "SE CCC" "Z" "SE Z"}; 

%if &type=d1 %then %do;
*** RHO1 ***;
PSI1_I = J(N,1,0);
PSI2_I = J(N,1,0);
DO I = 1 TO N;
   PSI1 = 0;
   PSI2 = 0;
   DO J = 1 TO N;

  %if &delta ne 0 %then %do;
   ** SUM OVER J **;
   B = 0.5*((abs(X[J,1:p]-Y[J,1:p])##&delta)*D1*(abs(X[J,1:p]-Y[J,1:p])##&delta)`);
   A = 0.5*((abs(X[I,1:p]-Y[J,1:p])##&delta)*D1*(abs(X[I,1:p]-Y[J,1:p])##&delta)`) + 0.5*((abs(X[J,1:p]-Y[I,1:p])##&delta)*D1*(abs(X[J,1:p]-Y[I,1:p])##&delta)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D1*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
   AI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D1*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
%END;


  %if &delta eq 0 %then %do;
	** SUM OVER J **;
	DIF11 = ((X[J,1:p]-Y[J,1:p])^=0);
	DIF12 = ((X[I,1:p]-Y[J,1:p])^=0);
	DIF13 = ((X[J,1:p]-Y[I,1:p])^=0);
	DIF14 = ((X[I,1:p]-Y[I,1:p])^=0);

   B = 0.5*((DIF11)*D1*(DIF11)`);
   A = 0.5*((DIF12)*D1*(DIF12)`) + 0.5*((DIF13)*D1*(DIF13)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (DIF14)*D1*(DIF14)`;
   AI = (DIF14)*D1*(DIF14)`;
	%END;

 

   *******************************;
   ** SUBTRACT OUT I=J **;
   PSI1_I[I] = PSI1 + ((N-2)/2)#BI;
   PSI2_I[I] = PSI2 - AI;

END;
U = SUM(PSI1_I)/(N*(N-1));
V = SUM(PSI2_I)/(N*(N-1));
** CCC IN TERMS OF U-STATISTICS **;
num1 = (n-1)*(v-u);
denom1 = u+(n-1)*v;
*print num1[label='Numerator'] denom1[label='Denominator'];
RHOCRM1 = ((N-1)*(V-U))/(U+(N-1)*V);
** CODE FOR ASYMPTOTIC VARIANCE **;
PSI1I = PSI1_I/(N-1);
PSI2I = PSI2_I/(N-1);
PSI_I = (PSI1I||PSI2I)`;
Uvec = U//V;
U_M = REPEAT(Uvec,1,N);
DELTA = PSI_I - U_M;
c = {1 2, 2 4};
V_U = (c/(N)##2)#(DELTA*DELTA`);
partialU = (-n#(n-1)#V)/(U+(n-1)#V)##2;
partialV = (n#(n-1)#U)/(U+(n-1)#V)##2;
d = partialU||partialV;
V_RHO = d*V_U*d`;
SE_RHO = SQRT(V_rho);
T=PROBIT(0.975);
Z=LOG((1+rhocrm1)/(1-rhocrm1))/2;
SE_Z = SQRT(V_rho)/(1-rhocrm1##2);
Z_LCL=Z-(SE_Z#T);
Z_UCL=Z+(SE_Z#T);
LCL1=(EXP(2#Z_LCL)-1)/(EXP(2#Z_LCL)+1);
UCL1=(EXP(2#Z_UCL)-1)/(EXP(2#Z_UCL)+1);

output=rhocrm1 || lcl1 || ucl1 || se_rho || Z || se_z; 
print output[label='' 
 			colname=labels rowname=' ' format=15.5];
print &type[label='Weight Matrix'];

%end; 

%else %if &type=d2 %then %do;
*** RHO2 ***;
PSI1_I = J(N,1,0);
PSI2_I = J(N,1,0);
DO I = 1 TO N;
   PSI1 = 0;
   PSI2 = 0;
   DO J = 1 TO N;
    %if &delta ne 0 %then %do;
   ** SUM OVER J **;
   B = 0.5*((abs(X[J,1:p]-Y[J,1:p])##&delta)*D2*(abs(X[J,1:p]-Y[J,1:p])##&delta)`);
   A = 0.5*((abs(X[I,1:p]-Y[J,1:p])##&delta)*D2*(abs(X[I,1:p]-Y[J,1:p])##&delta)`) + 0.5*((abs(X[J,1:p]-Y[I,1:p])##&delta)*D2*(abs(X[J,1:p]-Y[I,1:p])##&delta)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D2*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
   AI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D2*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
%END;


  %if &delta eq 0 %then %do;
	** SUM OVER J **;
	DIF11 = ((X[J,1:p]-Y[J,1:p])^=0);
	DIF12 = ((X[I,1:p]-Y[J,1:p])^=0);
	DIF13 = ((X[J,1:p]-Y[I,1:p])^=0);
	DIF14 = ((X[I,1:p]-Y[I,1:p])^=0);

   B = 0.5*((DIF11)*D2*(DIF11)`);
   A = 0.5*((DIF12)*D2*(DIF12)`) + 0.5*((DIF13)*D2*(DIF13)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (DIF14)*D2*(DIF14)`;
   AI = (DIF14)*D2*(DIF14)`;
	%END;
   *******************************;
   ** SUBTRACT OUT I=J **;
   PSI1_I[I] = PSI1 + ((N-2)/2)#BI;
   PSI2_I[I] = PSI2 - AI;
END;
U = SUM(PSI1_I)/(N*(N-1));
V = SUM(PSI2_I)/(N*(N-1));
** CCC IN TERMS OF U-STATISTICS **;
num2 = (n-1)*(v-u);
denom2 = u+(n-1)*v;
*print num2 denom2;
RHOCRM2 = ((N-1)*(V-U))/(U+(N-1)*V);
** CODE FOR ASYMPTOTIC VARIANCE **;
PSI1I = PSI1_I/(N-1);
PSI2I = PSI2_I/(N-1);
PSI_I = (PSI1I||PSI2I)`;
Uvec = U//V;
U_M = REPEAT(Uvec,1,N);
DELTA = PSI_I - U_M;
c = {1 2, 2 4};
V_U = (c/(N)##2)#(DELTA*DELTA`);
partialU = (-n#(n-1)#V)/(U+(n-1)#V)##2;
partialV = (n#(n-1)#U)/(U+(n-1)#V)##2;
d = partialU||partialV;
V_RHO = d*V_U*d`;
SE_RHO = SQRT(V_rho);
T=PROBIT(0.975);
Z=LOG((1+rhocrm2)/(1-rhocrm2))/2;
SE_Z = SQRT(V_rho)/(1-rhocrm2##2);
Z_LCL=Z-(SE_Z#T);
Z_UCL=Z+(SE_Z#T);
LCL2=(EXP(2#Z_LCL)-1)/(EXP(2#Z_LCL)+1);
UCL2=(EXP(2#Z_UCL)-1)/(EXP(2#Z_UCL)+1);

output=rhocrm2 || lcl2 || ucl2 || se_rho || Z || se_z; 
print output[label='' 
 			colname=labels rowname=' ' format=15.5];
print &type[label='Weight Matrix'];
%end;

%else %if &type=d3 %then %do;
*** RHO3 ***;
PSI1_I = J(N,1,0);
PSI2_I = J(N,1,0);
DO I = 1 TO N;
   PSI1 = 0;
   PSI2 = 0;
   DO J = 1 TO N;
    %if &delta ne 0 %then %do;
   ** SUM OVER J **;
   B = 0.5*((abs(X[J,1:p]-Y[J,1:p])##&delta)*D3*(abs(X[J,1:p]-Y[J,1:p])##&delta)`);
   A = 0.5*((abs(X[I,1:p]-Y[J,1:p])##&delta)*D3*(abs(X[I,1:p]-Y[J,1:p])##&delta)`) + 0.5*((abs(X[J,1:p]-Y[I,1:p])##&delta)*D3*(abs(X[J,1:p]-Y[I,1:p])##&delta)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D3*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
   AI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D3*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
%END;


  %if &delta eq 0 %then %do;
	** SUM OVER J **;
	DIF11 = ((X[J,1:p]-Y[J,1:p])^=0);
	DIF12 = ((X[I,1:p]-Y[J,1:p])^=0);
	DIF13 = ((X[J,1:p]-Y[I,1:p])^=0);
	DIF14 = ((X[I,1:p]-Y[I,1:p])^=0);

   B = 0.5*((DIF11)*D3*(DIF11)`);
   A = 0.5*((DIF12)*D3*(DIF12)`) + 0.5*((DIF13)*D3*(DIF13)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (DIF14)*D3*(DIF14)`;
   AI = (DIF14)*D3*(DIF14)`;
	%END;
   *******************************;
   ** SUBTRACT OUT I=J **;
   PSI1_I[I] = PSI1 + ((N-2)/2)#BI;
   PSI2_I[I] = PSI2 - AI;
END;
U = SUM(PSI1_I)/(N*(N-1));
V = SUM(PSI2_I)/(N*(N-1));
** CCC IN TERMS OF U-STATISTICS **;
num3 = (n-1)*(v-u);
denom3 = u+(n-1)*v;
RHOCRM3 = ((N-1)*(V-U))/(U+(N-1)*V);
** CODE FOR ASYMPTOTIC VARIANCE **;
PSI1I = PSI1_I/(N-1);
PSI2I = PSI2_I/(N-1);
PSI_I = (PSI1I||PSI2I)`;
Uvec = U//V;
U_M = REPEAT(Uvec,1,N);
DELTA = PSI_I - U_M;
c = {1 2, 2 4};
V_U = (c/(N)##2)#(DELTA*DELTA`);
partialU = (-n#(n-1)#V)/(U+(n-1)#V)##2;
partialV = (n#(n-1)#U)/(U+(n-1)#V)##2;
d = partialU||partialV;
V_RHO = d*V_U*d`;
SE_RHO = SQRT(V_rho);
T=PROBIT(0.975);
Z=LOG((1+rhocrm3)/(1-rhocrm3))/2;
SE_Z = SQRT(V_rho)/(1-rhocrm3##2);
Z_LCL=Z-(SE_Z#T);
Z_UCL=Z+(SE_Z#T);
LCL3=(EXP(2#Z_LCL)-1)/(EXP(2#Z_LCL)+1);
UCL3=(EXP(2#Z_UCL)-1)/(EXP(2#Z_UCL)+1);

output=rhocrm3 || lcl3 || ucl3 || se_rho || Z || se_z; 
print output[label='' 
 			colname=labels rowname=' ' format=15.5];
print &type[label='Weight Matrix'];
%end;

%else %if &type=d4 %then %do;
*** RHO4 ***;
PSI1_I = J(N,1,0);
PSI2_I = J(N,1,0);
DO I = 1 TO N;
   PSI1 = 0;
   PSI2 = 0;
   DO J = 1 TO N;
   %if &delta ne 0 %then %do;
   ** SUM OVER J **;
   B = 0.5*((abs(X[J,1:p]-Y[J,1:p])##&delta)*D4*(abs(X[J,1:p]-Y[J,1:p])##&delta)`);
   A = 0.5*((abs(X[I,1:p]-Y[J,1:p])##&delta)*D4*(abs(X[I,1:p]-Y[J,1:p])##&delta)`) + 0.5*((abs(X[J,1:p]-Y[I,1:p])##&delta)*D4*(abs(X[J,1:p]-Y[I,1:p])##&delta)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D4*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
   AI = (abs(X[I,1:p]-Y[I,1:p])##&delta)*D4*(abs(X[I,1:p]-Y[I,1:p])##&delta)`;
%END;


  %if &delta eq 0 %then %do;
	** SUM OVER J **;
	DIF11 = ((X[J,1:p]-Y[J,1:p])^=0);
	DIF12 = ((X[I,1:p]-Y[J,1:p])^=0);
	DIF13 = ((X[J,1:p]-Y[I,1:p])^=0);
	DIF14 = ((X[I,1:p]-Y[I,1:p])^=0);

   B = 0.5*((DIF11)*D4*(DIF11)`);
   A = 0.5*((DIF12)*D4*(DIF12)`) + 0.5*((DIF13)*D4*(DIF13)`);
   *******************************;
   PSI1 = PSI1 + B;
   PSI2 = PSI2 + A;
   END;
   ** DEFINE FUNCTIONS FOR I=J **;
   BI = (DIF14)*D4*(DIF14)`;
   AI = (DIF14)*D4*(DIF14)`;
	%END;
   *******************************;
   ** SUBTRACT OUT I=J **;
   PSI1_I[I] = PSI1 + ((N-2)/2)#BI;
   PSI2_I[I] = PSI2 - AI;
END;
U = SUM(PSI1_I)/(N*(N-1));
V = SUM(PSI2_I)/(N*(N-1));
** CCC IN TERMS OF U-STATISTICS **;
num4 = (n-1)*(v-u);
denom4 = u+(n-1)*v;
RHOCRM4 = ((N-1)*(V-U))/(U+(N-1)*V);
** CODE FOR ASYMPTOTIC VARIANCE **;
PSI1I = PSI1_I/(N-1);
PSI2I = PSI2_I/(N-1);
PSI_I = (PSI1I||PSI2I)`;
Uvec = U//V;
U_M = REPEAT(Uvec,1,N);
DELTA = PSI_I - U_M;
c = {1 2, 2 4};
V_U = (c/(N)##2)#(DELTA*DELTA`);
partialU = (-n#(n-1)#V)/(U+(n-1)#V)##2;
partialV = (n#(n-1)#U)/(U+(n-1)#V)##2;
d = partialU||partialV;
V_RHO = d*V_U*d`;
SE_RHO = SQRT(V_rho);
T=PROBIT(0.975);
Z=LOG((1+rhocrm4)/(1-rhocrm4))/2;
SE_Z = SQRT(V_rho)/(1-rhocrm4##2);
Z_LCL=Z-(SE_Z#T);
Z_UCL=Z+(SE_Z#T);
LCL4=(EXP(2#Z_LCL)-1)/(EXP(2#Z_LCL)+1);
UCL4=(EXP(2#Z_UCL)-1)/(EXP(2#Z_UCL)+1);

output=rhocrm4 || lcl4 || ucl4 || se_rho || Z || se_z; 
print output[label='' 
 			colname=labels rowname=' ' format=15.5];
print &type[label='Weight Matrix'];
%end;

create &output from output[colname=labels];
append from output;


%MEND UCCC;


****************************************************************************
CCC_LON macro produces the estimated CCC-RM via variance components for an 
identity or diagonal weight matrix. 
****************************************************************************;
%macro ccc_lon(dades,output,rho=0,D=0,detail=off);
%put NOTE: ccc_lon detail is: &detail.;
%if &detail=off or &detail=OFF or &detail=Off %then %do; ods select none; %end;



* Models with different options of within-subject correlation;

	proc mixed data=&dades asycov method=reml;
	class ind met time;
	model y = met time met*time/s covb;
	random int met time /subject=ind type=vc;
    %if (&rho=1) %then %do;	repeated time/subject=ind(met) type=ar(1) r; %end; /*Autoregressive*/
	%else %if (&rho=2) %then %do; repeated time/subject=ind(met) type=cs r; %end; /*Compound Symmetry*/
	lsmeans met*time /e;
	ods output SolutionF=_table1;
	ods output CovB=_table2;
	ods output CovParms=_table3;
	ods output AsyCov=_table4;
	ods output ClassLevels=_table5;
	ods output coef=_table6;
	run;

ods select all;
proc iml;
	use _table1; read all var{'Estimate'} into beta; 							
	use _table2; read all into aux;												
	use _table3; read all var('Estimate') into se where(CovParm='Residual');		
	use _table3; read all var('Estimate') into sa where(CovParm='Intercept');
	use _table3; read all var('Estimate') into sag where(CovParm='time');
	use _table3; read all var('Estimate') into sab where(CovParm='met');
	use _table4; read all into Saux;
	use _table5; read all var('Levels') into levels;
	use _table6; read all into auxC;

	ns=levels[1]; nm=levels[2]; nt=levels[3];  * ns,nm,nt stand for #subjects, #methods & #times;
	S1=Saux[,2:ncol(Saux)];
	
* var-covar fixed effects matrix;
	VarB=aux[,(1+ncol(aux)-nrow(aux)):ncol(aux)];

*building L matrix: design matrix of fixed effects;
C=auxC[,(1+ncol(auxC)-nt*nm):ncol(auxC)];


nd=nm*(nm-1)/2;					* Number of differences;
L=j(nrow(beta),nt*nd,0);
k=0;
do i=1 to nt*(nm-1);
	do j=1 to (nm-1);
		if ((i+nt*j)<=(nt*nm)) then do;
			k=k+1;
			L[,k]=C[,i]-C[,i+nt*j];
		end;
	end;
end;



**********************;
* Unweighted option **;
**********************;

%if (&D=0) %then %do;
* calculating sb: between-observers variability;

	difmed=t(L)*beta;							* vector of differences of means;
	A=L*t(L);

	aux1=(t(difmed)*difmed)-trace(A*VarB);	
	sb=max(aux1/(nm*(nm-1)*nt),0);	

* calculating the CCC;
	den=sa+sab+sag+se+sb;	
	ccc=(sa+sag)/den;							

* Variance of between-observers variability;

	var_sb=((2*trace((A*varB)**2))+(4*(t(beta)*A*varB*A*beta)))/((nm*(nm-1)*nt)**2); * Variance of sb;

*Covariance between sb and the remaining parameters;	

* dev: Vector of derivatives;

	dev_sa=(1-ccc)/den;
	dev_sag=(1-ccc)/den;
	dev_sb=(-1)*ccc/den;	
	if (sb=0) then dev_sb=0;
	dev_sab=(-1)*ccc/den;
	dev_se=(-1)*ccc/den;
	dev_rho=0;


	if (&rho=0) then do;
		cov_sasb=(-1/ns)*(S1[1,2]+S1[1,4]);	
		cov_sabsb=(-1/ns)*(S1[2,2]+S1[2,4]);	
		cov_sagsb=(-1/ns)*(S1[3,2]+S1[3,4]);	
		cov_sbrho=0;	
		cov_sbse=(-1/ns)*(S1[4,2]+S1[4,4]);

		S2 = cov_sasb || cov_sabsb || cov_sagsb || cov_sbse ||  var_sb;

		S=j(5,5,0);
		S[1:4,1:4]=S1;
		S[5,]=S2;
		S[1:4,5]=S2[1:4];

	end;
	else do;
		cov_sasb=(-1/ns)*(S1[1,2]+S1[1,5]);	
		cov_sabsb=(-1/ns)*(S1[2,2]+S1[2,5]);	
		cov_sagsb=(-1/ns)*(S1[3,2]+S1[3,5]);	
		cov_sbrho=0;	
		cov_sbse=(-1/ns)*(S1[5,2]+S1[5,5]);

		S2 = cov_sasb || cov_sabsb || cov_sagsb || cov_sbrho || cov_sbse ||  var_sb;

		S=j(6,6,0);
		S[1:5,1:5]=S1;
		S[6,]=S2;
		S[1:5,6]=S2[1:5];

	end;
%end;


**********************;
* Weighted option   **;
**********************;


%if (&D^=0) %then %do;
	use &D; read all into D;
	auxD=D;
		if (nd > 1) then do;
				auxD=block(D,D);
				c=2;
			do while (c<nd);
				c=c+1;
				auxD=block(auxD,D);
			end;
	
		end;
	

		difmed=t(L)*beta;							* vector of differences of means;
		AW=L*auxD*t(L);
		aux1=(t(difmed)*auxD*difmed)-trace(AW*VarB);	
			
		sb=max(aux1/(nm*(nm-1)),0);
		sumd=sum(D);
		
		* calculating the CCC;

		den=(sumd*(sa+sab+sag+se))+sb;	
		ccc=(sumd*(sa+sag))/den;					*weighted ccc;	


	* Variance of between-observers variability;
		var_sb=((2*trace((AW*varB)**2))+(4*(t(beta)*AW*varB*AW*beta)))/((nm*(nm-1))**2);


		dev_sa=sumd*(1-ccc)/den;
		dev_sag=sumd*(1-ccc)/den;
		dev_sb=(-1)*ccc/den;	
		if (sb=0) then dev_sb=0;
		dev_sab=sumd*(-1)*ccc/den;
		dev_rho=0;
		dev_se=sumd*(-1)*ccc/den;


		if (&rho=0) then do;

			var_se=S1[4,4];							* variance of random error variance;
			aux2=trace(AW*(varB/se));

			cov_sasb=(-1/ns)*(S1[1,2]+S1[1,4]);	
			cov_sabsb=(-1/ns)*(S1[2,2]+S1[2,4]);	
			cov_sagsb=(-1/ns)*(S1[3,2]+S1[3,4]);	
			cov_sbrho=0;	
			cov_sbse=(-1/ns)*(S1[4,2]+S1[4,4]);

			S2 = cov_sasb || cov_sabsb || cov_sagsb || cov_sbse ||  var_sb;

			S=j(5,5,0);
			S[1:4,1:4]=S1;
			S[5,]=S2;
			S[1:4,5]=S2[1:4];

		end;
		else do;

			var_se=S1[5,5];							* variance of random error variance;
			aux2=trace(AW*(varB/se));

			cov_sasb=(-1/ns)*(S1[1,2]+S1[1,5]);	
			cov_sabsb=(-1/ns)*(S1[2,2]+S1[2,5]);	
			cov_sagsb=(-1/ns)*(S1[3,2]+S1[3,5]);	
			cov_sbrho=0;	
			cov_sbse=(-1/ns)*(S1[5,2]+S1[5,5]);

			S2 = cov_sasb || cov_sabsb || cov_sagsb || cov_sbrho || cov_sbse ||  var_sb;

			S=j(6,6,0);
			S[1:5,1:5]=S1;
			S[6,]=S2;
			S[1:5,6]=S2[1:5];

		end;
%end;


%if (&rho=0) %then %do;
	dev= dev_sa || dev_sab || dev_sag ||  dev_se || dev_sb;
%end;
%else %do;
	dev= dev_sa || dev_sab || dev_sag ||  dev_rho || dev_se || dev_sb;	
%end;

* Confidence interval estimation;
	%ic_ccc(ccc,dev,S);
	
* Output;
	varcomp= sa || sab || sag || sb || se;

	print 'Repeated Measures Concordance Correlation Coefficient','Derived from Variance Components',,
        ccc[label='' colname=('CCC')];
	print 'Variance Components',,
		varcomp[colname={'Subjects' 'Subjects-Observer' 'Subjects-Time' 'Observers' 'Error'} 
		label='' format=20.4 rowname=''];
	print '','','CCC Output',,result[colname={'CCC' 'Lower 95% CL' 'Upper 95% CL' 'SE CCC' 'Z' 'SE Z'} rowname='' label='' format=15.4];
	%if (&D^=0) %then %do; print D[label='Weight Matrix']; %end;

	create &output from result[colname={'VC_CCC' 'VC_LL95' 'VC_UL95' 'VC_SECCC' 'VC_Z' 'VC_SEZ'}];
	append from result;

%exit:

quit;

%mend;


***********************************************************************************************
CCC estimates the CCC through variance components from a linear mixed model with 
subject and observer effects (no subject by observer interaction) and observers considered as 
a fixed effect if m= the number of measures by subject and observer is equal to 1.  Otherwise, 
the macro estimates the CCC through variance components including the subject by observer 
interaction;
*************************************************************************************************;
%macro ccc(infile,output=ccc,m=1,detail=off);

%put NOTE: ccc detail is: &detail.;
%if &detail=off or &detail=OFF or &detail=Off %then %do; ods select none; %end;

proc mixed data=&infile asycov method=reml;
	class met ind;
	model y=met/covb solution;
	%if &m = 1 %then %do; random ind ; %end;
    %else %if &m > 1 %then %do; random ind ind*met ; %end;
	lsmeans met ;
	ods output SolutionF=_table1;
	ods output CovB=_table2;
	ods output CovParms=_table3;
	ods output AsyCov=_table4;
	ods output ClassLevels=_table5;
	run;

ods select all;
proc iml;
	use _table1; read all var{'Estimate'} into betaaux; 							
	use _table2; read all into aux;
	use _table3; read all var('Estimate') into VC;		
	use _table4; read all into Saux;
	use _table5; read all var('Levels') into levels;
	S1=Saux[,2:ncol(Saux)];
	sa=VC[1]; 
    %if &m=1 %then %do; se=VC[2]; %end;
    %else %if &m>1 %then %do; sab=VC[2]; se=VC[3]; %end;
	n=levels[2]; k=levels[1]; ncat=nrow(levels);
	
* var-covar fixed effects matrix;
    VarB=aux[,ncol(aux)-k:ncol(aux)]; 

*building L matrix;

	* aux1 = design matrix;	
	aux1=j(k+1,k,0);
	aux1[1,1:k]=1;
		
	do i=1 to (k-1);	aux1[1+i,1+i]=1; end;
		
	L=j(k+1,k*(k-1)/2,0);
	cont=0;	
	do i=1 to k;
		do j=1 to (k-1) while (i>j);
			cont=cont+1;
			L[,cont]=aux1[,i]-aux1[,j];
		end;
	end;
			
* calculating sb;
	beta=betaaux[1:(k+1)]; 
	difmed=t(L)*beta;							* vector of differences of means;
	var_difmed=t(L)*varB*L;
	A=L*t(L);
	aux2=(t(difmed)*difmed)-trace(A*VarB);	
	sb=max(aux2/(k*(k-1)),0);	

* CCC;
	den=sa+se+sb;	
	ccc=sa/den;							

* Standard Error;

	S1=Saux[,2:ncol(Saux)];
	var_sb=((2*trace((A*varB)**2))+(4*(t(beta)*A*varB*A*beta)))/(&m*(k*(k-1))**2);
	
	%if &m>1 %then %do;	
	cov_sasb=(-1/n*&m)*S1[1,3];	
	cov_sabsb=(-1/n*&m)*S1[2,3];	
	cov_sbse=(-1/n*&m)*S1[3,3];
	
	S2 = cov_sasb || cov_sabsb || cov_sbse ||  var_sb;
	
	S=j(4,4,0);
	S[1:3,1:3]=S1;
	S[4,]=S2;
	S[1:4,4]=S2[1:4];

	* Derivatives;
	dev_sa=(1-ccc)/den;
	dev_sb=(-1)*ccc/den;	
	%if (sb=0) %then %do; dev_sb=0; %end;
	dev_sab=(-1)*ccc/den;
	dev_se=(-1)*ccc/den;
	
	dev= dev_sa  || dev_sab || dev_se || dev_sb;

	%end;

	%else %if &m=1 %then %do;
	cov_sasb=(-1/n*&m)*S1[1,2];	
	cov_sbse=(-1/n*&m)*S1[2,2];
	
	S2 = cov_sasb || cov_sbse ||  var_sb;
	
	S=j(3,3,0);
	S[1:2,1:2]=S1;
	S[3,]=S2;
	S[1:3,3]=S2[1:3];

	* Derivatives;

	dev_sa=(1-ccc)/den;
	dev_sb=(-1)*ccc/den;	
	%if (sb=0) %then %do; dev_sb=0; %end;
	dev_se=(-1)*ccc/den;
	
	dev= dev_sa  || dev_se || dev_sb;
    %end;

	%ic_ccc(ccc,dev,S);
	
	reset spaces = 3;
	%if &m=1 %then %do; varcomp= sa || se || sb; %end;
	%else %if &m>1 %then %do; varcomp = sa || sab || se || sb; %end;
	print 'Repeated Measures Concordance Correlation Coefficient','Derived from Variance Components',,
        ccc[label='' colname=('CCC')];
	%if &m>1 %then %do; print 'Variance Components',,varcomp[colname={'Subject' 'Subject by Observer' 'Random Error' 'Observer' } 
		format=20.8 label='' rowname='']; %end;
	%else %if &m=1 %then %do; print 'Variance Components',,varcomp[colname={'Subject' 'Random Error' 'Observer' } 
		format=20.8 label='' rowname='']; %end;
	print S[label="Variance-Covariance Matrix of Variance Components" format=20.8];
	res=t(result[1:6]);
	print 'Concordance Correlation Coefficient',,res[colname={'CCC' 'Lower 95% CL' 'Upper 95% CL' 'Standard Error CCC'
		  'Z' 'SE Z'} format=6.4 
	    rowname='' label=''];

	
	create &output from result[colname={'VC_CCC' 'VC_LL95' 'VC_UL95' 'VC_SECCC' 'VC_Z' 'VC_SEZ'}];
	append from result;
	quit;
%mend;

****************************************************************************;
* Confidence interval for CCC using Z tranformation;
* CCC is the estimate;
* dev is the vector of derivatives of the variance components to the CCC;
* S is the variance-covariance matrix of the variance components;
****************************************************************************;
%macro ic_ccc(CCC,dev,S);
	SECCC = sqrt(&dev*&S*&dev`);
	Z=0.5*log((1+&CCC)/(1-&CCC));
	SEZ=sqrt((SECCC**2)/(((1+&CCC)**2)*((1-&CCC)**2)));
	T=PROBIT(0.975);
	LLZ=Z-T*SEZ;
	ULZ=Z+T*SEZ;
	LL95 = (exp(2*LLZ)-1)/(exp(2*LLZ)+1);
	UL95 = (exp(2*ULZ)-1)/(exp(2*ULZ)+1);
	result= &CCC || LL95 || UL95 || SECCC || Z || SEZ;
%mend;


