function [loops,v] = PSOR_solver(v,g,b,Nminus,Nplus,c1,c2,c3,omega,eps)

% length of v
N = Nplus - Nminus + 1;
loops = 0;
error = 1;

while error >= eps
    
    error = 0.0;
    
    for n = 2:N-1    
    y = (b(n)-c1*v(n-1)-c3*v(n+1))/c2;
    y = max(g(n),v(n)+omega*(y-v(n)));
    error = error + (v(n)-y)^2;
    v(n) = y;   
    end
    
    loops = loops + 1;
    
end 

end
