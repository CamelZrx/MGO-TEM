function z = store(x_ini,xo,eo,L)
% x_ini is the extremum seed, xo is the design vector, eo is the function
% value of the design vector and L is the level information of the extremum
% seed. z is the newly stored local minimum
global ObSet EAP Su_MP W1t cs ct Etrue CurrentBest pop
global num Ebest nonlcon_i tol_eq_i 
global num_consec_best prvs_best num_consec_best_all fmin fmax
persistent ee 

if(isempty(ee))
    ee=1;
end

if(norm(x_ini-xo.*cs-ct)~=0)    
    ObSet(length(ObSet)+1).seed=x_ini;      % store into the observation set
    ObSet(length(ObSet)).vector=xo;
    ObSet(length(ObSet)).value=eo;
    ObSet(length(ObSet)).level=L;
    ObSet(length(ObSet)).radius=norm(x_ini-xo.*cs-ct);    
end

if(~isempty(nonlcon_i))
    [c,ceq] = nonlcon_i(ObSet(end).vector);
else
    c = 0;
    ceq = 0;
end

if(CurrentBest.value>eo && sum(c>0)==0 && sum(abs(ceq)>tol_eq_i)==0)     % check feasibility
    CurrentBest=ObSet(end);
end

if(CurrentBest.value<prvs_best)
    num_consec_best = 0;
    prvs_best = CurrentBest.value;
else
    num_consec_best = num_consec_best + 1;    
    prvs_best = CurrentBest.value;
end
num_consec_best_all = [num_consec_best_all num_consec_best];

pop(end).value = CurrentBest.value;
pop(end).vector = CurrentBest.vector;
Ebest(num) = CurrentBest.value;
Etrue(num) = eo;

if(eo<fmin)
    fmin = eo;                            % update the negmin if smaller negative local minimum is found
end

if(eo>fmax)
    fmax = eo;                            % update the negmin if smaller negative local minimum is found
end

[Np] = normalize_up_to_t(eo);

EAP = [EAP (eo-2*min(0,fmin)+1)/Np];                        % store into the normalized EAP set    
if(length(EAP)>=2)    
    w = exp(-min(EAP));
    Su_MP = [Su_MP w];                        % store into the sub-martingale set for nonnegative local minimums    
    if(length(Su_MP)<=1)
        W1t = [W1t 1];                        
    else
        W1t = [W1t 1-mean(Su_MP)/Su_MP(end-1)];         
    end
end

z = ObSet(end);

end