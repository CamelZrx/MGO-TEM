function z = tunneling(x1,x2,L)
% x1 is the first sample point and x2 is the second sample point.
global ObSet dim stopflg epsilon delta x_ini_tun cs ct Etrue maxover CurrentBest

N = length(ObSet);

if(stopflg(end)==1)
    z = x1;
    return
end

% expected too high
if(mean([x1.value,x2.value])>mean(Etrue))
    z = CurrentBest;
    return;
end

if(x1.value<=x2.value)
    Range = coordinate_cons(x1.vector.*cs+ct,x2.vector.*cs+ct);% create the subdomain    
else    
    Range = coordinate_cons(x2.vector.*cs+ct,x1.vector.*cs+ct); % create the subdomain   
end

for i =1:dim                            % ensure the sampling is within the box
    if(Range(i).ub>1)
        Range(i).ub=1;
    end
    if(Range(i).lb<0)
        Range(i).ub=0;
    end
end

for i =1:dim
    xt(i) = rand*(Range(i).ub-Range(i).lb)+Range(i).lb;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% may avoid solving if consider convex hull
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xo,eo,maxover] = local_m(xt);          % solve the optimization problem using the 
z = store(xt,xo,eo,L+1);                % tunneling extremum seed and store it
sfg = check_stop(eo,epsilon,delta);
x_ini_tun = [x_ini_tun;xt L+1];

stopflg = [stopflg sfg];

if(sfg==1)
    z = x1;
    return
end

if(maxover==1)
    z = x1;
    return;
end

end