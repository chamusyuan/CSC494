function [loops,x] = Policy_Iterator(A,x0,b,c,Nminus,Nplus)

loops = 0;
error = 1.0;
%expected tolerance
eps = 10^-8;

N = Nplus - Nminus - 1;
%x is the old v(m+1)
x = x0;

phi = zeros(N,1);
oldPHI = ones(N,1);
newPHI = zeros(N,1);

%stop the iteration when error is negligible or phi doesn't change
while (error > eps && (isequal(oldPHI,newPHI)==0))

    %disp(loops)
    oldPHI = newPHI;
    
    newA = sparse(N,N);
    newb = zeros(N,1);
    
    % find optimal phi
    for i = 1:N
    phi0 = x(i) - c(i);
    temp = A*x;
    phi1 = temp(i) - b(i);
    
    if phi1 < phi0
        phi(i) = 1;
    end
    
    %update A* and b*
    identity = speye(N); 
    temp = sparse(i,1:N,phi(i)*A(i,:).' + (1-phi(i))*identity(i,:).',N,N);
    newA = newA + temp;
    newb(i) = phi(i)*b(i) + (1-phi(i))*c(i);
    end
    
    newPHI = phi;

    %solve new v(m+1)
    newx = newA \ newb;
    
    error = max(abs(newx - x));
    x = newx;
    loops = loops + 1;

end

end