
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

%compute front and rear cornerning stifness
Ca_r=(a/L)*m*g*B*C*D;
Ca_f=(b/L)*m*g*B*C*D;
% 
%2.1.1 compute the front and rear cornering stifness for the vehicle generate equilibrium trajetory using Euler integration and linear tire
%forces

delta_fun=@(t) 10*pi/180*sin(2*pi*t)-0.0175;
Z_eq = zeros(5,101);
Z_eq(:,1) = [0 0 0 0 0]';
x(1) = Z_eq(1,1);
y(1) = Z_eq(2,1);
psi(1) = Z_eq(3,1);
vy(1) = Z_eq(4,1);
r(1) = Z_eq(5,1);

for i = 1:1:length(T)-1
    alpha_f = delta_fun(T(i)) - (vy(i)+a*r(i))/vx;
    alpha_r = -(vy(i)-b*r(i))/vx;
    Fyf = Ca_f*alpha_f;
    Fyr = Ca_r*alpha_r;

    xdot = vx*cos(psi(i)) - vy(i)*sin(psi(i));
    ydot = vy(i)*cos(psi(i)) + vx*sin(psi(i));
    psidot = r(i);
    vydot = (1/m)*(Fyr + Fyf - m*vx*r(i));
    rdot = (1/Iz)*(-b*Fyr + a*Fyf);
    
    x(i+1) = x(i) + (T(i+1) - T(i))*xdot;
    y(i+1) = y(i) + (T(i+1) - T(i))*ydot;
    psi(i+1) = psi(i) + (T(i+1) - T(i))*psidot;
    vy(i+1) = vy(i) + (T(i+1) - T(i))*vydot;
    r(i+1) = r(i) + (T(i+1) - T(i))*rdot;
end

Z_eq=[x; y; psi; vy; r];


%2.1.2 linearization for feedback gains bike with linear tire forces
AFun = @(i)[0, 0, -vx*sin(psi(i))-vy(i)*cos(psi(i)), -sin(psi(i)), 0;
     0, 0, -vy(i)*sin(psi(i))+vx*cos(psi(i)), cos(psi(i)), 0;
     0, 0, 0, 0, 1;
     0, 0, 0, (-Ca_r-Ca_f)/(vx*m), (Ca_r*b-Ca_f*a)/(vx*m)-vx;
     0, 0, 0, (Ca_r*b-Ca_f*a)/(vx*Iz), (-Ca_r*b^2-Ca_f*a^2)/(vx*Iz)];

BFun = @(i)[0;
     0;
     0;
     Ca_f/m;
     a*Ca_f/Iz];
 
Q = eye(5);
R = 0.5;
[K, P] = lqr_LTV(AFun,BFun,Q,R,T);

%2.1.3 Plot linear vs nonlinear tire forces and find max % difference

for i = 1:1:length(T)
    alphaFmagic(i) = delta_fun(T(i))-atan((vy(i)+a*r(i))/vx);
    alphaRmagic(i) = -atan((vy(i)-b*r(i))/vx);
    FyfMagic(i) = (b*m*g/L)*D*sin(C*atan(B*(1-E)*alphaFmagic(i)+E*atan(B*alphaFmagic(i))));
    FyrMagic(i) = (a*m*g/L)*D*sin(C*atan(B*(1-E)*alphaRmagic(i)+E*atan(B*alphaRmagic(i))));
    alpha_f(i) = delta_fun(T(i)) - (vy(i)+a*r(i))/vx;
    alpha_r(i) = -(vy(i)-b*r(i))/vx;
    Fyf(i) = Ca_f*alpha_f(i);
    Fyr(i) = Ca_r*alpha_r(i);
end

FrontWheelError = 100*max(abs(Fyf - FyfMagic)./abs(FyfMagic));
RearWheelError = 100*max(abs(Fyr - FyrMagic)./abs(FyrMagic));

tireforce_percent_error = max(FrontWheelError, RearWheelError);

%2.1.4 Euler Simulate with  Nonlinear tire dynamics

Z = zeros(5,101);
Z(:,1) = [0 0 0 0 0]';
deltaLin = @(Z,i) K{i}*(Z_eq(:,i) - Z(:,i)) + delta_fun(T(i));
xL(1) = Z_eq(1,1);
yL(1) = Z_eq(2,1);
psiL(1) = Z_eq(3,1);
vyL(1) = Z_eq(4,1);
rL(1) = Z_eq(5,1);

for i = 1:1:length(T)-1
    deltaL = deltaLin(Z, i);
    if deltaL > 45*pi/180
        deltaL = 45*pi/180;
    elseif deltaL < -45*pi/180
        deltaL = -45*pi/180;
    end
    alphaFmagic(i) = deltaL-atan((vyL(i)+a*rL(i))/vx);
    alphaRmagic(i) = -atan((vyL(i)-b*rL(i))/vx);
    FyfMagic(i) = (b*m*g/L)*D*sin(C*atan(B*(1-E)*alphaFmagic(i)+E*atan(B*alphaFmagic(i))));
    FyrMagic(i) = (a*m*g/L)*D*sin(C*atan(B*(1-E)*alphaRmagic(i)+E*atan(B*alphaRmagic(i))));
    
    xdot = vx*cos(psiL(i)) - vyL(i)*sin(psiL(i));
    ydot = vyL(i)*cos(psiL(i)) + vx*sin(psiL(i));
    psidot = rL(i);
    vydot = (1/m)*(FyrMagic(i) + FyfMagic(i)*cos(deltaL) - m*vx*rL(i));
    rdot = (1/Iz)*(-b*FyrMagic(i) + a*FyfMagic(i)*cos(deltaL));
    
    xL(i+1) = xL(i) + (T(i+1) - T(i))*xdot;
    yL(i+1) = yL(i) + (T(i+1) - T(i))*ydot;
    psiL(i+1) = psiL(i) + (T(i+1) - T(i))*psidot;
    vyL(i+1) = vyL(i) + (T(i+1) - T(i))*vydot;
    rL(i+1) = rL(i) + (T(i+1) - T(i))*rdot;
    Z(:,i+1) = [xL(i+1); yL(i+1); psiL(i+1); vyL(i+1); rL(i+1)];
end    

max_distance_error=max(sqrt((x - xL).^2 + (y - yL).^2));

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