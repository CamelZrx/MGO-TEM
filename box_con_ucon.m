function [z1,z2,z3,z4] = box_con_ucon(lb,ub)
% lb is the input vector for the lower bounds, ub is the input vector for
% the upper bounds; z1 is the contraction scale factor, z2 is the
% contraction translation; z3 is the uncontraction scale factor, z4 is the
% uncontraction translation.
global dim

    for i = 1:dim
        scale(i) = 1/(ub(i)-lb(i));       % contraction scale factor
        trans(i) = -lb(i)/(ub(i)-lb(i));  % contraction translation
    end
    
    z1 = scale;
    z2 = trans;
    z3 = 1./scale;                                  % uncontraction scale
    z4 = -trans.*z3;                                % uncontraction translation
end