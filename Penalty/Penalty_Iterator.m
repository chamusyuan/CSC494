function [value] = Penalty_Iterator(oldv,g,b,A,tol,N,L)

    loops = 0;
    % initial setting for penalty matrix
    oldP = sparse(1:N,1:N,L*(oldv<g),N,N);
    
    %nonlinear euqation f
    f = A*oldv + b + oldP*(g-oldv);
    % first order derivative of f at v(m+1)
    df = sparse(A - oldP);
    
    % solving for new v(m+1)
    newit = oldv - df \ f;
    oldit = oldv;
    
    while max(abs(newit - oldit)) > tol
        %disp(loops)
        oldit = newit;
        oldP = sparse(1:N,1:N,L*(oldit<g),N,N);
        
        %updating f
        f = A*oldit + b  + oldP*(g-oldit);
        %updating df
        df = sparse(A - oldP);
        %updating v(m+1)
        newit = oldit - df \ f;
        %updating penalty matrix
        newP = sparse(1:N,1:N,L*(newit<g),N,N);
        %stop the iteration if penalty doesn't change
        if newP == oldP
            break
        end
        loops = loops+1;
    end
    value = newit;
end