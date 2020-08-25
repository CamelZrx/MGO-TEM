function [z1,z2,z3]=local_m(x_ini)
% x_ini is the extremum seed, z1 is the design vector, z2 is the function
% value and z3 is the maximum overflow flag
global dim ucs uct Numfeval MAX opts objfun_i A_i b_i Aeq_i beq_i nonlcon_i

x_ini_real = x_ini.*ucs+uct;            %transfer back to the original box
lb_real = zeros(1,dim).*ucs+uct;
ub_real = ones(1,dim).*ucs+uct;

[x,fval,~,output] = fmincon(objfun_i,x_ini_real,A_i,b_i,Aeq_i,beq_i,lb_real,ub_real,nonlcon_i,opts);

Numfeval = [Numfeval Numfeval(end)+output.funcCount];

z1 = x;                           % design vector
z2 = fval;

nn = Numfeval(end);
if(nn>MAX)
    z3=1;
else
    z3=0;
end

end