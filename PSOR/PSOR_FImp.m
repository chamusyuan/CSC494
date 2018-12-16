function [values] = PSOR_FImp(dx,dt,M,Nplus,Nminus,volatility,r)
%Implicit FD & PSOR for American Option

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

%c1,c2,c3 are coefficients w.r.t v(n-1,m+1),v(n,m+1),v(n+1,m+1)
c1 = - 0.5*dt*volatility^2/dx^2 + 0.5*dt*(r-0.5*volatility^2)/dx;
c2 = dt*volatility^2/dx^2 + dt*r + 1;
c3 = - 0.5*dt*volatility^2/dx^2 - 0.5*dt*(r-0.5*volatility^2)/dx;

b = zeros(N,1);
for j = 1:M
    
    %b is the right hand side of the system to solve
    for n = 2:(N-1)
        b(n) = v(n);
    end
    
    %boundary conditions
    v(1) = 1;
    v(N) = 0;
    
    [loops,v] = PSOR_solver(v,g,b,Nminus,Nplus,c1,c2,c3,omega,eps);
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


