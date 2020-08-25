function z = coordinate_cons(x1,x2)
% x1 is the first extremum seed; x2 is the second extremum seed.

global dim;

z = struct('lb',[],'ub',[]);

distt = norm(x1-x2,2);       % compute the distance between x1 and x2
len = distt/sqrt(dim);       % length of the hyper cubic
c = (x1+x2)/2;              % center of the hyper ball
for i = 1:dim
    z(i).lb = max(0,c(i)-len/2);
    z(i).ub = min(1,c(i)+len/2);
end

end