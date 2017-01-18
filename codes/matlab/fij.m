function [f,gx,gy] = fij(alpha,tau,lij,uij,x,y)
%FIJ
%    OUT1 = FIJ(ALPHA,LIJ,TAU,UIJ,X1,X2,X3,Y1,Y2,Y3)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    17-Jan-2017 13:34:12

x1 = x(1); x2 = x(2); x3 = x(3);
y1 = y(1); y2 = y(2); y3 = y(3); 
t2 = x1-y1;
t3 = x2-y2;
t4 = x3-y3;
t6 = tau.^2;
t7 = t2.^2;
t8 = t3.^2;
t9 = t4.^2;
t10 = t6+t7+t8+t9;
t11 = sqrt(t10);
t5 = -t11+uij;
t12 = alpha.^2;
t13 = lij-t11;
t14 = t13.^2;
t15 = t12.*t14;
t16 = t6+t15;
t17 = x1.*2.0;
t29 = y1.*2.0;
t18 = t17-t29;
t19 = t5.^2;
t20 = t12.*t19;
t21 = t6+t20;
t22 = 1.0./sqrt(t10);
t23 = 1.0./sqrt(t16);
t24 = x2.*2.0;
t30 = y2.*2.0;
t25 = t24-t30;
t26 = 1.0./sqrt(t21);
t27 = x3.*2.0;
t31 = y3.*2.0;
t28 = t27-t31;
f  =  -alpha.*t5+alpha.*t13+sqrt(t16)+sqrt(t21);
gx = [t5.*t12.*t18.*t22.*t26.*(-1.0./2.0)-t12.*t13.*t18.*t22.*t23.*(1.0./2.0);
      t5.*t12.*t22.*t25.*t26.*(-1.0./2.0)-t12.*t13.*t22.*t23.*t25.*(1.0./2.0);
      t5.*t12.*t22.*t26.*t28.*(-1.0./2.0)-t12.*t13.*t22.*t23.*t28.*(1.0./2.0)];
gy = [t5.*t12.*t18.*t22.*t26.*(1.0./2.0)+t12.*t13.*t18.*t22.*t23.*(1.0./2.0);
      t5.*t12.*t22.*t25.*t26.*(1.0./2.0)+t12.*t13.*t22.*t23.*t25.*(1.0./2.0);
      t5.*t12.*t22.*t26.*t28.*(1.0./2.0)+t12.*t13.*t22.*t23.*t28.*(1.0./2.0)];
