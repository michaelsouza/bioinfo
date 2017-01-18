function out1 = fij_grad_x(alpha,lij,tau,uij,x1,x2,x3,y1,y2,y3)
%FIJ_GRAD_X
%    OUT1 = FIJ_GRAD_X(ALPHA,LIJ,TAU,UIJ,X1,X2,X3,Y1,Y2,Y3)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    17-Jan-2017 13:29:00

t2 = alpha.^2;
t3 = x1-y1;
t4 = x2-y2;
t5 = x3-y3;
t7 = tau.^2;
t8 = t3.^2;
t9 = t4.^2;
t10 = t5.^2;
t11 = t7+t8+t9+t10;
t12 = sqrt(t11);
t6 = lij-t12;
t13 = x1.*2.0;
t14 = t13-y1.*2.0;
t15 = t12-uij;
t16 = 1.0./sqrt(t11);
t17 = t6.^2;
t18 = t2.*t17;
t19 = t7+t18;
t20 = 1.0./sqrt(t19);
t21 = x2.*2.0;
t22 = t21-y2.*2.0;
t23 = t12-uij;
t24 = x3.*2.0;
t25 = t24-y3.*2.0;
t26 = t12-uij;
out1 = [t2.*t14.*t16.*1.0./sqrt(t7+t2.*t15.^2).*(t12-uij).*(1.0./2.0)-t2.*t6.*t14.*t16.*t20.*(1.0./2.0);t2.*t16.*t22.*1.0./sqrt(t7+t2.*t23.^2).*(t12-uij).*(1.0./2.0)-t2.*t6.*t16.*t20.*t22.*(1.0./2.0);t2.*t16.*t25.*1.0./sqrt(t7+t2.*t26.^2).*(t12-uij).*(1.0./2.0)-t2.*t6.*t16.*t20.*t25.*(1.0./2.0)];
