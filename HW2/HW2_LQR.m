d=1.5;
l=3;
gam=atan(0.3);
v=2;
psi=@(t) (v/l*tan(gam))*t;

% 1.1 Simulate the non-linear system to generate the circular trajectory
f=@(t,x) [v/l*(l*cos(psi(t))-d*sin(psi(t))*tan(gam));
          v/l*(l*sin(psi(t))+d*cos(psi(t))*tan(gam));
          v/l*tan(gam)];
T= 0:.01:5;

[T,Y]=ode45(f,T,[0,0,0]);

%1.2 Find the linearized time varying A and B matricies
A=@(t) [0 0 (v/l)*(-l*sin(.2*t)-d*cos(.2*t)*tan(gam));
        0 0 (v/l)*(l*cos(.2*t)-d*sin(.2*t)*tan(gam));
        0 0 0];
    
B=@(t) [1/l*(l*cos(.2*t)-d*sin(.2*t)*tan(gam)) -d*(v/l)*sin(.2*t)/cos(gam)^2;
        1/l*(l*sin(.2*t)+d*cos(.2*t)*tan(gam)) d*(v/l)*cos(.2*t)/cos(gam)^2;
        1/l*tan(gam) v/l*1/cos(gam)^2];

%1.3 Find the optimal feedback gains
Q=eye(3);
R=eye(2);

%the function lqr_LTV requires the A and B matrix functions to be functions of the step number i
A1=@(i) A(T(i));
B1=@(i) B(T(i));
[K,P]=lqr_LTV(A1,B1,Q,R,T);

%1.4
x0=[0.1;0.8;0.01];

%similarly to lqr_LTV this function requires f to be a function of the step i
f1=@(i,dx) A1(i)*dx - B1(i)*K{i}*dx;
Y1=ode1(f1,T,x0)+Y;