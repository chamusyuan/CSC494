function [values] = Penalty_Method_CN(dx,dt,M,Nplus,Nminus,volatility,r,tol,L)
%crank-nicolson & Penalty Method for American Option

N = Nplus - Nminus - 1;

%Setting initial values for v and g, notice that they stand at tau = 0, 
%which means t = T, and g = oldv.
oldv = zeros(N,1);
g = zeros(N,1);
for i = 1:N
    oldv(i) = max(1-exp((Nminus+i)*dx),0);
    g(i) = oldv(i);
end

%z4,z5,z6 are coefficients w.r.t v(n-1,m),v(n,m),v(n+1,m)
z1 = dt * volatility^2 / dx^2 - dt*(r-0.5*volatility^2) / dx;
z2 = -4 - 2*dt*volatility^2 / dx^2 - 2*dt*r;
z3 = dt * volatility^2 / dx^2 + dt*(r-0.5*volatility^2) / dx;

%z1,z2,z3 are coefficients w.r.t v(n-1,m+1),v(n,m+1),v(n+1,m+1)
z4 = z1;
z5 = 4 - 2*dt*volatility^2 / dx^2 - 2*dt*r;
z6 = z3;

%constructing left hand side matrix A
D = sparse(1:N,1:N,z2,N,N);
E = sparse(2:N,1:N-1,z1,N,N);
F = sparse(1:N-1,2:N,z3,N,N);
A = E+D+F;

%constructing right hand side matrix B
G = sparse(1:N,1:N,z5,N,N);
H = sparse(2:N,1:N-1,z4,N,N);
I = sparse(1:N-1,2:N,z6,N,N);
B = G+H+I;

%boundary conditions
first = sparse(1,1,1,N,1);
last = sparse(N,1);

for j = 1:M
 
    %b is the right hand side of the system to solve
    b = z1*first + z3*last + B*oldv + z4*first + z6*last;    
    
    newv = Penalty_Iterator(oldv,g,b,A,tol,N,L);
    oldv = newv;
    
for i = 1:N
    values(i) = oldv(i);
end
    
end