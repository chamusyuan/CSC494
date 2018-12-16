function [values] = PSOR_CN(dx,dt,M,Nplus,Nminus,volatility,r)
%crank-nicolson & PSOR for American Option

N = Nplus - Nminus + 1;
%expected tolerance
eps = 10^(-16);
omega = 1.0;
domega = 0.05;
oldloops = 10000;

%Setting initial values for v and g, notice that they stand at tau = 0, 
%which means t = T, and g = v.
v = zeros(N,1);
g = zeros(N,1);
for i = 1:N
    v(i) = max(1-exp((Nminus+i-1)*dx),0);
    g(i) = v(i);
end

%z4,z5,z6 are coefficients w.r.t v(n-1,m),v(n,m),v(n+1,m)
z4 = dt * volatility^2 / dx^2 - dt*(r-0.5*volatility^2) / dx;
z5 = 4 - 2*dt*volatility^2 / dx^2 - 2*dt*r;
z6 = dt * volatility^2 / dx^2 + dt*(r-0.5*volatility^2) / dx;

%z1,z2,z3 are coefficients w.r.t v(n-1,m+1),v(n,m+1),v(n+1,m+1)
z1 = -z4;
z2 = 4 + 2*dt*volatility^2 / dx^2 + 2*dt*r;
z3 = -z6;

b = zeros(N,1);

for j = 1:M
    
    %b is the right hand side of the system to solve
    for n = 2:(N-1)
        b(n) = z4*v(n-1) + z5*v(n) + z6*v(n+1);
    end
    
    %boundary conditions
    v(1) = 1;
    v(N) = 0; 
    
    [loops,v] = PSOR_solver(v,g,b,Nminus,Nplus,z1,z2,z3,omega,eps);
    %changing omega
    if (loops > oldloops)
        domega = domega*(-1.0);
    end
    omega = omega + domega;
    oldloops = loops;
    
end

%output values    
for i = 1:N
    values(i) = v(i);
end
    
end


