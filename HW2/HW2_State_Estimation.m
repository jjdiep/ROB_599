%parameter definition
ms=400;
A=[0 -1;0 0];
B=[0;-1/ms];
F=[1;0];
C=[1 0];

%1.1 determine the observability matrix (O) and its rank (r)
O= [C; C*A];
r= 2;
%1.2 Pole placement Design
%1.2.1 Feedback Design. Find gain matrix, K, to place poles at -1,-1 
K=[0 1]*inv([B A*B])*(A*A+2*A+eye(2));
%1.2.2 Observer Design. Find gain matrix, G, to place poles at -3,-3 
G=(A^2+6*A+9*eye(2))*O'*[0; 1];

%1.2.3 Euler simulate to find the trajectories for 
%the system (x) and observer (x_hat)
t=0:1/100:10;
x_1=[0.1;0.1];
x_hat_1=[0;0];

x = zeros(2,length(t));
x_hat = zeros(2,length(t));

x(:,1) = x_1;
x_hat(:,1) = x_hat_1;

for i = 1:length(t)-1
    ode_x = A*x(:,i)-B*K*x_hat(:,i)+[randn(1)*sqrt(.0005);0];
    ode_xhat = (A-G*C)*x_hat(:,i)+G*(C*x(:,i)+randn(1)*sqrt(.0005))-B*K*x_hat(:,i);
    x(:,i+1) = x(:,i) + (t(i+1) - t(i))*ode_x;
    x_hat(:,i+1) = x_hat(:,i) + (t(i+1) - t(i))*ode_xhat;
end

x = x;
x = x_hat;


%1.3 LQR Design
%1.3.1 Feedback Design. use LQR to find feedback gains K_lqr
%note to use the lqr_LTV function with a constant A, B matrices, 
%enter @(i)A, @(i)B as the arguments AFun, BFun
Q=100*eye(2);
R=0.00005;

K_lqr= lqr_LTV(@(i)A,@(i)B,Q,R,t);

%1.3.2 Observer Design. use LQR to find optimal observer G_lqr
%note to use the lqr_LTV function with a constant A, B matrices, 
%enter @(i)A, @(i)B as the arguments AFun, BFun
Qo=eye(2);
Ro=1;
[G_trans,P_2] = lqr_LTV(@(i)A',@(i)C',Qo,Ro,t);
G_lqr= G_trans;
G_lqr = cellfun(@transpose,G_lqr,'UniformOutput',false);

%1.3.3 Euler Simulate

x(:,1) = x_1;
x_hat(:,1) = x_hat_1;

for i = 1:length(t)-1
    ode_x = A*x(:,i)-B*cell2mat(K_lqr(i))*x_hat(:,i)+[randn(1)*sqrt(.0005);0];
    ode_xhat = (A-cell2mat(G_lqr(i))*C)*x_hat(:,i)+cell2mat(G_lqr(i))*(C*x(:,i)+randn(1)*sqrt(.0005))-B*cell2mat(K_lqr(i))*x_hat(:,i);
    x(:,i+1) = x(:,i) + (t(i+1) - t(i))*ode_x;
    x_hat(:,i+1) = x_hat(:,i) + (t(i+1) - t(i))*ode_xhat;
end

x_lqr= x;
x_hat_lqr= x_hat;

%1.4 Comparison
%Look at the observer error (x-x_hat) and (x_lqr-x_hat_lqr)


%which observer converges faster? enter 'acker' or 'lqr'
%apostrophes included!
faster_convergence='acker';

%which observer is less noisy? enter 'acker' or 'lqr'
%apostrophes included!
less_noisy='lqr';

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