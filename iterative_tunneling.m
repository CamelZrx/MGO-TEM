function z = iterative_tunneling(x1,x2,L)
% x1 is the firtst sample point and x2 is the second sample point, L is the
% level information of the input sample point. z is the tunneling point at
% the Lth level.
global cs tol1 tol2 ObSet dim CurrentBest stopflg

if(stopflg(end)==1)
    z = ObSet(end);
    return
end

if(L>sqrt(dim))    
    z = ObSet(end);
    return
end

xt = tunneling(x1,x2,L);                        % get the tunneling extremum seed
if(norm(xt.vector.*cs-x1.vector.*cs)/sqrt(dim)<tol1&&abs(xt.value-x1.value)<tol2)
                                                % no new local minimum in the subdomain
elseif(norm(xt.vector.*cs-x2.vector.*cs)/sqrt(dim)<tol1&&abs(xt.value-x2.value)<tol2) 
                                                % no new local minimum in the subdomain
else
    if(xt.value==CurrentBest.value)
%         disp('exploitation')
        if(x1.value<x2.value)
            iterative_tunneling(CurrentBest,x1,L+1);    % further tunneling required
        else
            iterative_tunneling(CurrentBest,x2,L+1);    % further tunneling required
        end        
    else
%         disp('tunnel')
        iterative_tunneling(x1,xt,L+1);    % further tunneling required        
    end
end

z = ObSet(end);

end