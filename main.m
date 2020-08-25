clc
options = optimset('LargeScale','off','MaxFunEvals',500*10,...
    'MaxIter',1000,'Algorithm','active-set','Display','off');

tic
dim = 5;
box1 = -10;
box2 = 10;
maxeval = inf;500*dim;
epsilon = 0.2;
delta = 0.06;
tol = 1e-3;
[a,b] = MDGOP(@(x) objfun(x),dim,[],[],[],[],box1*ones(dim,1),box2*ones(dim,1),...
    [],maxeval,epsilon,delta,options,tol,'Plot');
toc
a
b(end)