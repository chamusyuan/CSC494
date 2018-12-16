function [values] = Penalty_Method_FImp(dx,dt,M,Nplus,Nminus,volatility,r,tol,L)
%fully-implicit & Penalty Method for American Option

N = Nplus - Nminus - 1;

%Setting initial values for v and g, notice that they stand at tau = 0, 
%which means t = T, and g = oldv.
oldv = zeros(N,1);
g = zeros(N,1);
for i = 1:N
    oldv(i) = max(1-exp((Nminus+i)*dx),0);
    g(i) = oldv(i);
end

%c1,c2,c3 are coefficients w.r.t v(n-1,m+1),v(n,m+1),v(n+1,m+1)
c1 = 0.5*dt*volatility^2/dx^2 - 0.5*dt*(r-0.5*volatility^2)/dx;
c2 = -dt*volatility^2/dx^2 - dt*r - 1;
c3 = 0.5*dt*volatility^2/dx^2 + 0.5*dt*(r-0.5*volatility^2)/dx;

%constructing left hand side matrix A
D = sparse(1:N,1:N,c2,N,N);
E = sparse(2:N,1:N-1,c1,N,N);
F = sparse(1:N-1,2:N,c3,N,N);
A = E+D+F;

b = zeros(N,1);

%boundary conditions
first = sparse(1,1,1,N,1);
last = sparse(N,1);

for j = 1:M

    %b is the right hand side of the system to solve
    b = oldv + c1*first + c3*last;
    newv = Penalty_Iterator(oldv,g,b,A,tol,N,L);
    oldv = newv;
    
for i = 1:N
    values(i) = oldv(i);
end
    
end


