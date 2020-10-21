%%
A=[-3 0 0;0 -1 -3;0 1 -2];
[eigvect,eigval] = eig(A);
%2.1
l1= eigval(1,1);
l2= eigval(2,2);
l3= eigval(3,3);

%List eigenvectors such that v1 is the eigenvector for eigenvalue l1.
%Eigenvectors should be in column vector format.
v1=[eigvect(:,1)];
v2=[eigvect(:,2)];
v3=[eigvect(:,3)];

%2.2
%Write out the close form equation using the symbolic variable t as time
syms t
x_0=[5;0;0];
x=expm(A*t)*x_0;

%2.3
f=@(t,x) A*x;
tspan=[0:.1:10];
x0=[0; 1; 2];

[T,Y]=ode45(f,tspan,x0);

%2.4 
sol24=0;