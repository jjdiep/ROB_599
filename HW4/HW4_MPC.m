% Problem 1.1
A_cont = [0 .310625; 0 0];
B_cont = [0 ; 1];
dt = .01;
A = eye(2) + A_cont*dt;
B = dt*B_cont;

% Problem 1.2
Ndec = 11+11+10;

% Problem 1.3
Aeq = zeros(22,32);
beq = zeros(22,1);
beq(1) = -1.25-(-1);

x = zeros(2,11);
x(:,1) = [1; 0];
Aeq(1:2,1:2) = eye(2);
for k = 3:2:21
    Aeq(k:k+1,k-2:k+1) = [A -eye(2)];
end
for k = 3:2:21
    Aeq(k:k+1,(k+1)/2+21) = B;
end
% Problem 1.4
u = @(t) -0.3175*sin(pi*t/10-pi/2);

Aineq = [eye(32); -eye(32)];
bineq = [ones(22,1)*.5;ones(10,1)*(10-u(0));ones(22,1)*.5;ones(10,1)*(10+u(0))];
% Problem 1.5
H = [eye(22)*100 zeros(22,10); zeros(10,32)];
c = zeros(32,1);
dx_k = [-.25; 0];
x_k = [-1.25; 0];

beq2 = beq;
% bineq2 = bineq;

for i = 1:1:1000-1
    solDV = quadprog( H, c, Aineq, bineq, Aeq, beq2 );
    du_k = solDV(23);
    u_nom = u(.01*(i-1));
    
    dx_k(:,i+1) = A*dx_k(:,i)+B*du_k;
    x_k(:,i+1) = A*x_k(:,i)+B*(du_k+u_nom);
    
    beq2(1:2) = dx_k(:,i+1);
%     bineq2 = [ones(22,1)*.5;ones(10,1)*(10-u(.01*0));ones(22,1)*.5;ones(10,1)*(10+u(.01*0))];
end
x = x_k';
% transpose x at the end
