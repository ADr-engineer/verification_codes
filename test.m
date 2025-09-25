clc
clear
% 
% Ef=7e4;
% Es=7e4;
% vf=0.17;
% vs=0.17;
% hf=3;
% hs=10;
% Cf=(1-vf^2)/(Ef*hf);
% Cs=-(1-vs^2)/(Es*hs);
% af=3.8e-5;
% as=3.5e-5;
% R=50;
% T1=500;
% T2=-5e-2;
% T3=-2.5e-1;
% r=-50:0.1:50;

Ef=7e10;
Es=7e10;
vf=0.17;
vs=0.17;
hf=3e-3;
hs=10e-3;
Cf=(1-vf^2)/(Ef*hf);
Cs=-(1-vs^2)/(Es*hs);
af=3.8e-5;
as=3.5e-5;
R=0.05;
T1=500;
T2=-5e4;
T3=-2.5e5;
r=0:0.001:0.05;


An = [-Cs/(4*(Cf-Cs)*hf), 3*Cs/(2*(Cf-Cs)*hf^2), (4*Cf-3*Cs)/(4*(Cf-Cs)*hs), 3*Cs/(2*(Cf-Cs)*hs^2);
      -3*Cf/(2*(Cf-Cs)*hf^2), -3*(Cf-4*Cs)/((Cf-Cs)*hf^3), 3*Cf/(2*(Cf-Cs)*hf*hs), 9*Cf/((Cf-Cs)*hf*hs^2);
      3*Cs/(2*(Cf-Cs)*hf*hs), -9*Cs/((Cf-Cs)*hf^2*hs), -3*Cs/(2*(Cf-Cs)*hs^2), -3*(4*Cf-Cs)/((Cf-Cs)*hs^3)];

Yr_1n = [(1+vf)/Cf -(1+vf)*hf/(2*Cf) -hs*(1+vf)/(2*Cf)];
Yr_2n = [(3*Cf-4*Cs)*(-1+vf)/(4*Cf*(Cf-Cs)*hf),3*(-1+vf)/(2*(Cf-Cs)*hf^2),...
         (-1+vf)/(4*(Cf-Cs)*hs), 3*(-1+vf)/(2*(Cf-Cs)*hs^2)];
Yr_3n = [-1/(4*(Cf-Cs)*hf), 3/(2*(Cf-Cs)*hf^2), 3/(2*(Cf-Cs)*hs),...
         3/(2*(Cf-Cs)*hs^2)];

Yr_1u = [(1+vf)/(Cf*hf) -(1+vf)/Cf -hs*(1+vf)/(2*Cf*hf)];
Yr_2u = [(3*Cf-2*Cs)*(-1+vf)/(2*Cf*(Cf-Cs)*hf^2),3*(Cf-2*Cs)*(-1+vf)/(Cf*(Cf-Cs)*hf^3),...
         (1-vf)/(2*(Cf-Cs)*hf*hs), -3*(-1+vf)/((Cf-Cs)*hf*hs^2)];
Yr_3u = [(3*Cf-2*Cs)/(2*Cf*(Cf-Cs)*hf^2), 3*(Cf-2*Cs)/(Cf*(Cf-Cs)*hf^3),...
         -1/(2*(Cf-Cs)*hf*hs), -3/((Cf-Cs)*hf*hs^2)];
Ytheta_1u = Yr_1u;
Ytheta_2u = -Yr_2u;
Ytheta_3u = vf*Yr_3u;

Yr_1m = [(1+vf)/(Cf*hf) -(1+vf)/(2*Cf) -hs*(1+vf)/(2*Cf*hf)];
Yr_2m = Yr_2n/hf;
Yr_3m = Yr_2m/(vf-1);
Ytheta_1m = Yr_1m;
Ytheta_2m = -Yr_2m;
Ytheta_3m = vf/(1-vf)*Ytheta_2m;


Yr_1l = [(1+vf)/(Cf*hf) 0 -hs*(1+vf)/(2*(Cf*hf))];
Yr_2l = [Cs*(1-vf)/(Cf*(Cf-Cs)*hf^2),6*Cs*(-1+vf)/(Cf*(Cf-Cs)*hf^3),...
         (-1+vf)/((Cf-Cs)*hf*hs), 6*(-1+vf)/((Cf-Cs)*hf*hs^2)];
Yr_3l = Yr_2l/(vf-1);
Ytheta_1l = Yr_1l;
Ytheta_2l = -Yr_2l;
Ytheta_3l = vf/(1-vf)*Ytheta_2l;


b11 = Cs*(Cs*(1-vf)+Cf*(-1+vs))/(4*(Cf-Cs)*hf*R^2*(-Cs*(1+vf)+Cf*(1+vs)));
b12 = 3*Cs*(Cs*(-1+vf)+Cf*(1-vs))/(2*(Cf-Cs)*hf^2*R^2*(-Cs*(1+vf)+Cf*(1+vs)));
b13 = (-3*Cs^2*(1+vf)*(-1+vs)-4*Cf^2*(-1+vs^2)+Cf*Cs*(-7+2*vs+3*vs^2+vf*(-2+4*vs))) / (4*(Cf-Cs)*hs*R^2*(1+vs)*(-Cs*(1+vf)+Cf*(1+vs)));
b14 = 3*Cs*(Cs*(1+vf)*(-1+vs)+Cf*(1+2*vf-2*vs-vs^2))/(2*(Cf-Cs)*hs^2*R^2*(1+vs)*(-Cs*(1+vf)+Cf*(1+vs)));
b21 = 3*Cf*(-Cs*(-1+2*vf-2*vs+vf^2)+Cf*(-1+vf)*(1+vs))/(2*(Cf-Cs)*hf^2*R^2*(1+vf)*(-Cs*(1+vf)+Cf*(1+vs)));
b22 = (12*Cs^2*(-1+vf^2)+3*Cf^2*(-1+vf)*(1+vs)-3*Cf*Cs*(-5+2*vs+vf^2+vf*(-2+4*vs))) /((Cf-Cs)*hf^3*R^2*(1+vf)*(-Cs*(1+vf)+Cf*(1+vs)));
b23 = Cf*hf/(Cs*hs)*b12;
b24 = 6*Cf*hf/(Cs*hs^2)*b12;
b31 = hf/hs*b12;
b32 = -6/hs*b12;
b33 = -b14;
b34 = (3*Cs^2*(1+vf)*(-1+vs)+12*Cf^2*(-1+vs^2)-3*Cf*Cs*(-5-2*vs+vs^2+vf*(2+4*vs))) /((Cf-Cs)*hf^3*R^2*(1+vf)*(-Cs*(1+vf)+Cf*(1+vs)));

Bfree = [b11 b12 b13 b14;
         b21 b22 b23 b24;
         b31 b32 b33 b34];

z1 = hs/2+hf;
z2 = hs/2;
z3 = -hs/2;

A1 = (1+vf)*af*(T1*R^2/2*(z1-z2)+T2*R^4*(z1-z2)/4+T3*R^2*(z1^3-z2^3)/6);
A2 = (1+vf)*af*((1/4*T1*R^2*(z1^2-z2^2) + 1/8*T2*R^4*(z1^2-z2^2) + 1/8*T3*R^2*(z1^4-z2^4)) - ...
    (hs+hf)/2*(1/2*T1*(z1-z2)*R^2 + 1/4*T2*R^4*(z1-z2) + 1/6*T3*R^2*(z1^3-z2^3)));
A3 = (1+vs)*as*(1/2*T1*R^2*(z2-z3) + 1/4*T2*R^4*(z2-z3) + 1/6*T3*R^2*(z2^3-z3^3));
A4 = (1+vs)*as*(1/4*T1*R^2*(z2^2-z3^2) + 1/8*T2*R^4*(z2^2-z3^2) + 1/8*T3*R^2*(z2^4-z3^4));
Ac = [A1 A2 A3 A4]';



F1r = [(1+vf)*af*(T1/2*(z1-z2)+T2*r.^2*(z1-z2)/4+T3*(z1^3-z2^3)/6);
        (1+vf)*af*((1/4*T1*(z1^2-z2^2) + 1/8*T2*r.^2*(z1^2-z2^2) + 1/8*T3*(z1^4-z2^4)) - ...
    (hs+hf)/2*(1/2*T1*(z1-z2) + 1/4*T2*r.^2*(z1-z2) + 1/6*T3*(z1^3-z2^3)));
    (1+vs)*as*(1/2*T1*(z2-z3) + 1/4*T2*r.^2*(z2-z3) + 1/6*T3*(z2^3-z3^3));
    (1+vs)*as*(1/4*T1*(z2^2-z3^2) + 1/8*T2*r.^2*(z2^2-z3^2) + 1/8*T3*(z2^4-z3^4))];


F2r = [(1+vf)*af*(T1*(z1-z2)+T2*r.^2*(z1-z2)+1/3*T3*(z1^3-z2^3));
        (1+vf)*af*(1/2*T1*(z1^2-z2^2)+1/2*T2*r.^2*(z1^2-z2^2)+1/4*T3*(z1^4-z2^4)-(hs+hf)/2*...
        (T1*(z1-z2)+T2*r.^2*(z1-z2)+1/3*T3*(z1^3-z2^3)));
        (1+vs)*as*(T1*(z2-z3)+T2*r.^2*(z2-z3)+1/3*T3*(z2^3-z3^3));
        (1+vs)*as*(1/2*T1*(z2^2-z3^2)+1/2*T2*r.^2*(z2^2-z3^2)+1/4*T3*(z2^4-z3^4))];

st_r_up = Yr_1u*Bfree*Ac + Yr_2u*F1r + Yr_3u*F2r - (1+vf)*af*(T1+T2*r.^2+T3*z1^2)/(Cf*hf);
st_theta_up = Ytheta_1u*Bfree*Ac + Ytheta_2u*F1r + Ytheta_3u*F2r - (1+vf)*af*(T1+T2*r.^2+T3*z1^2)/(Cf*hf);
st_r_nt = Yr_1m*Bfree*Ac + Yr_2m*F1r + Yr_3m*F2r - (1+vf)*af*(T1+T2*r.^2+T3*(hs/2+hf/2)^2)/(Cf*hf);
st_theta_nt = Ytheta_1m*Bfree*Ac + Ytheta_2m*F1r + Ytheta_3m*F2r - (1+vf)*af*(T1+T2*r.^2+T3*(hs/2+hf/2)^2)/(Cf*hf);
st_r_bt = Yr_1l*Bfree*Ac + Yr_2l*F1r + Yr_3l*F2r - (1+vf)*af*(T1+T2*r.^2+T3*z2^2)/(Cf*hf);
st_theta_bt = Ytheta_1l*Bfree*Ac + Ytheta_2l*F1r + Ytheta_3l*F2r - (1+vf)*af*(T1+T2*r.^2+T3*z2^2)/(Cf*hf);


plot(r,st_r_nt,'o-')

save('variables.mat','An','F1r','F2r','Bfree','Yr_1n','Yr_2n','Yr_3n','Yr_1u',...
    'Yr_2u','Yr_3u','Yr_1m','Yr_1l','Yr_2l')
