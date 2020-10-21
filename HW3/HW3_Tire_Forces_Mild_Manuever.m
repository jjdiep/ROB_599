%Vehicle Parameters
a   =  1.14;		% distance c.g. to front axle (m) 
L   =  2.54;		% wheel base (m)
m   =  1500;		% mass (kg)
Iz  =  2420.0;	% yaw moment of inertia (kg-m^2)
b=L-a;   %distance of c.g to rear axel (m) 
g=9.81;
vx=20;

%%Tire forces
B=10;
C=1.3;
D=1;
E=0.97;


%timespan for all simulations
T=0:0.01:1;

%1.1 compute front and rear cornerning stifness
Ca_r=(a/L)*m*g*B*C*D;
Ca_f=(b/L)*m*g*B*C*D;
% 
%1.2.1 compute the front and rear cornering stifness for the vehicle generate equilibrium trajetory using Euler integration and linear tire
%forces

delta_fun=@(t) pi/180*sin(2*pi*t)-0.00175;
Z_eq = zeros(5,101);
Z_eq(:,1) = [0,0,0,0,0]';

x(1) = Z_eq(1,1);
y(1) = Z_eq(2,1);
psi(1) = Z_eq(3,1);
vy(1) = Z_eq(4,1);
r(1) = Z_eq(5,1);

for i = 1:1:length(T)-1
    alphaf = delta_fun(T(i))-(vy(i)+a*r(i))/vx;
    alphar = -(vy(i)-b*r(i))/vx;
    Fyf = Ca_f*alphaf;
    Fyr = Ca_r*alphar;
    
    xdot = vx*cos(psi(i))-vy(i)*sin(psi(i));
    ydot = vy(i)*cos(psi(i))+vx*sin(psi(i));
    psidot = r(i);
    vydot = 1/m*(Fyr+Fyf-m*vx*r(i));
    rdot = 1/Iz*(-b*Fyr+a*Fyf);
    
    x(i+1) = x(i)+(T(i+1)-T(i))*xdot;
    y(i+1) = y(i)+(T(i+1)-T(i))*ydot;
    psi(i+1) = psi(i)+(T(i+1)-T(i))*psidot;
    vy(i+1) = vy(i)+(T(i+1)-T(i))*vydot;
    r(i+1) = r(i)+(T(i+1)-T(i))*rdot;
             
end
Z_eq=[x;y;psi;vy;r];

%1.2.2 linearization for feedback gains bike with linear tire forces
A_fun = @(i) [0,0,-vx*sin(psi(i))-vy(i)*cos(psi(i)),-sin(psi(i)),0;
          0,0,-vy(i)*sin(psi(i))+vx*cos(psi(i)),cos(psi(i)),0;
          0,0,0,0,1;
          0,0,0,1/m*(-Ca_r/vx-Ca_f/vx),1/m*(Ca_r*(b/vx)+Ca_f*(-a/vx)-m*vx);
          0,0,0,1/Iz*(b*Ca_r/vx-a*Ca_f/vx),1/Iz*(-b*Ca_r*(b/vx)+a*Ca_f*(-a/vx))];
B_fun = @(i)[0;
         0;
         0;
         1/m*Ca_f;
         1/Iz*a*Ca_f];
Qo = eye(5);
R_val = 0.5;
[K,P] = lqr_LTV(A_fun,B_fun,Qo,R_val,T);

%1.2.3 Plot linear vs nonlinear tire forces and find max % difference

for i = 1:1:length(T)
    alphaf_mgc(i) = delta_fun(T(i))-atan((vy(i)+a*r(i))/vx);
    alphar_mgc(i) = -atan((vy(i)-b*r(i))/vx);
    Fyf_mgc(i) = (b/L*m*g)*D*sin(C*atan(B*(1-E)*alphaf_mgc(i)+E*atan(B*alphaf_mgc(i))));
    Fyr_mgc(i) = (a/L*m*g)*D*sin(C*atan(B*(1-E)*alphar_mgc(i)+E*atan(B*alphar_mgc(i))));
    
    alphaf_lin(i) = delta_fun(T(i))-(vy(i)+a*r(i))/vx;
    alphar_lin(i) = -(vy(i)-b*r(i))/vx;
    Fyf_lin(i) = Ca_f*alphaf_lin(i);
    Fyr_lin(i) = Ca_r*alphar_lin(i);
end
tireforce_percent_error_front = 100*max(abs(Fyf_lin-Fyf_mgc)./abs(Fyf_mgc));
tireforce_percent_error_rear = 100*max(abs(Fyr_lin-Fyr_mgc)./abs(Fyr_mgc));
tireforce_percent_error= max([tireforce_percent_error_front, tireforce_percent_error_rear]);

%1.2.4 Euler Simulate with  Nonlinear tire dynamics
Z = zeros(5,101);
Z(:,1) = [0,0,0,0,0]';

x_l(1) = Z(1,1);
y_l(1) = Z(2,1);
psi_l(1) = Z(3,1);
vy_l(1) = Z(4,1);
r_l(1) = Z(5,1);

delta = @(Z, i) K{i}*(Z_eq(:,i)-Z(:,i))+delta_fun(T(i));

for i = 1:1:length(T)-1
    delta_f = delta(Z, i);

    if delta_f > 45*pi/180
        delta_f = 45*pi/180;
    elseif delta_f < -45*pi/180
        delta_f = -45*pi/180;
    end
    
    alphaf = delta_f-atan((vy_l(i)+a*r_l(i))/vx);
    alphar = -atan((vy_l(i)-b*r(i))/vx);
    Fyf = (b/L*m*g)*D*sin(C*atan(B*(1-E)*alphaf+E*atan(B*alphaf)));
    Fyr = (a/L*m*g)*D*sin(C*atan(B*(1-E)*alphar+E*atan(B*alphar)));
    
    xdot = vx*cos(psi_l(i))-vy_l(i)*sin(psi_l(i));
    ydot = vy_l(i)*cos(psi_l(i))+vx*sin(psi_l(i));
    psidot = r_l(i);
    vydot = 1/m*(Fyr+Fyf*cos(delta_f)-m*vx*r_l(i));
    rdot = 1/Iz*(-b*Fyr+a*Fyf*cos(delta_f));
    
    x_l(i+1) = x_l(i)+(T(i+1)-T(i))*xdot;
    y_l(i+1) = y_l(i)+(T(i+1)-T(i))*ydot;
    psi_l(i+1) = psi_l(i)+(T(i+1)-T(i))*psidot;
    vy_l(i+1) = vy_l(i)+(T(i+1)-T(i))*vydot;
    r_l(i+1) = r_l(i)+(T(i+1)-T(i))*rdot;
    Z(:,i+1) = [x_l(i+1);y_l(i+1);psi_l(i+1);vy_l(i+1);r_l(i+1)];  
end
max_distance_error = max(sqrt((x-x_l).^2+(y-y_l).^2));
%Function library


function [K, P] = lqr_LTV(AFun,BFun,Q,R,tSpan)
    nSteps = length(tSpan);

    P{nSteps} = zeros(size(Q));
    K{nSteps} = zeros(length(R),length(Q));
    
    for i = nSteps-1:-1:1
        A_ = AFun(i+1);
        B_ = BFun(i+1);
        P_ = P{i+1};
        
        P{i} = P_ + (tSpan(i+1)-tSpan(i)) * ( P_*A_ + A_'*P_ - P_*B_*(R\(B_'*P_)) + Q);
        K{i} = R\(B_'*P_);
    end
end