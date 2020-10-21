%parameter definition
m=1500;
n=4.8;
r=0.4;
g=9.81;
C_r=0.01;
rho=1.3;
C_d=0.32;
a=2.4;
theta_e=2*pi/180;
v_e=20;

%you may enter equations, but do not change the names of the variables given. 
%The auto-grader uses these variable names to check your answers!

%Part 1.1 Determine equilibrium engine torque
u_e= (m*g*sign(v_e)*C_r+.5*rho*C_d*a*v_e^2+m*g*sin(theta_e))*r/n;

%Part 1.2.1 Determine the A,B,F matrices for the linearized system
A= -1/m*rho*C_d*a*v_e;
B= n/(m*r);
F= -g*cos(theta_e);

%Part 1.2.2 Determine feedback gains of linear system to place pole at -1
k= (A+1)/B;

%Part 1.2.3 Determine the steady state error after 10 seconds
theta_tilda=3*pi/180;
odefun = @(t,v_tilda) (A)*v_tilda+(B)*(-k)*v_tilda+F*theta_tilda;
v_0 = -1;
tspan = 0:.05:10;
[tout,v_tildaout] = ode45(odefun,tspan,v_0);
sse= v_tildaout(end);


%Part 1.3.1 Determine A_I,B_I,F_I matrices of the linear system with integral action
%the subscript I is just used to indicate these matrices and vectors apply 
A_I= [A 0; 1 0];
B_I= [B; 0];
F_I= [F; 0];

%Part 1.3.2 Place poles of the system with integral action at -1,-2
k_I= [0 1]*inv([B_I A_I*B_I])*(A_I*A_I+3*A_I+2*eye(2));

%Part 1.3.3 with integral action determine the steady state error after 10
%seconds
odefun2 = @(t2,v_line) [(A_I*v_line+B_I*(-k_I*v_line)+F_I*theta_tilda)];
v_0 = [-1; 0];
tspan2 = 0:.05:10;
[tout,v_lineout] = ode45(odefun2,tspan2,v_0);
sse_with_integral_action= v_lineout(1,end);
