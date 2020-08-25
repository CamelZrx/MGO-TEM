function [z1,z2] = MDGOP(fun,nvar,A,b,Aeq,beq,lb,ub,nonlcon,max_feval,eps,del,options,tol_eq,PlotFcn,ini_seed)
% fun is the objective function, nvar is the number of variables,lb are the lower bounds, ub are the upper
% bounds, max_feval is the maximum allowed function evaluation calls, eps
% is the epsilon confidence level for mean process and delta is the
% confidence level for sub-martingales, options are options specified for fmincon.
% ini_seed is the initial seed, specified by the user.
% z1 is the final best solution and z2 is the number of function evaluations used.

    global ObSet EAP Su_MP W1t epsilon delta Nt MP CMP Np
    global mt stopflg pflg nflg Numfeval dim cs ct ucs uct tol1 tol2 LB UB 
    global Etrue count x_inii x_ini_rand x_ini_tun x_ini_perm tol11 tol22 MAX
    global maxover CurrentBest BestSol pop opts objfun_i A_i b_i Aeq_i beq_i nonlcon_i
    global plotflag num maxrun Ebest UpP Up tol_eq_i num_consec_best prvs_best seed_i
    global num_consec_best_all alpha_tt beta_tt Etrue_T stop11 stop2a1 stop2b1 fmin fmax
    
    if(nargin<15)
        plotflag = 0;
    else
        if(strcmp(PlotFcn,'Plot'))
            plotflag = 1;
        else
            plotflag = 0;
        end
    end
    
    flag = 0;
    dim = nvar;                % dimension of the problem, i.e., the number of design variables
    num = 1;
    tol1 = 1e-3;
    tol2 = 1e-3;
    tol11 = 1e-3;
    tol22 = 1e-3;

    ObSet = struct('seed',[],'vector',[],'value',[],'level',[],'radius',[],'radiusfz',[]);
    ObSet(1)=[];
    CurrentBest = struct('seed',[],'vector',[],'value',[],'level',[],'radius',[],'radiusfz',[]);
    CurrentBest(1)=[];
    pop = struct('seed',[],'vector',[],'value',[],'level',[],'radius',[],'radiusfz',[]);
    pop = [];
    EAP = [];               % the experiment set for nonnegative local minimums
    Etrue = [];             % the true value of feval
    Ebest = [];             % the best value of feval    
    Su_MP = [];             % the sub-martingale of the nonnegative local minimums
    Nt = 0;                 % the normalization factor up to time t
    MP = [];                % the array for mean process
    UpP = 1;                % upcrossing for positive process
    Up = [];                % upcrossing history
    CMP = [];               % the cesaro mean process
    W1t = [];               % the stopping criteria two for the positive sub-martingale
    Np = 1;                 % normalization factor for nonnegative branch
    mt = 0;                 % statistics for sub-martingale history
    stopflg = 0;            % the stop flag, stop if raised to 1
    pflg = 0;               % the stop flag for nonnegative sub-martingale
    nflg = 0;               % the stop flag for negative sub-martingale
    Numfeval = 0;           % number of function evaluation calls
    count = 0;              % counting number for randomize
    epsilon = eps;          % 1-epsilon sure we have the true RV
    delta = del;           % 1-delta sure the sub-martingale is stable
    x_inii = rand(1,dim);   % test for sampling method
    x_ini_rand = [];        % Gibbs Sampling
    x_ini_tun = [];         % Gibbs Sampling in tunneling
    x_ini_perm = [];         % Gibbs Sampling in tunneling
    BestSol = [];
    MAX = max_feval;           % the maximum allowed function evaluation call
    maxover = 0;            % maximum function evaluation call over flag
    A_i = A;
    b_i = b;
    Aeq_i = Aeq;
    beq_i = beq;
    LB = lb';
    UB = ub';
    nonlcon_i = nonlcon;
    opts = options;
    objfun_i = fun;
    tol_eq_i = tol_eq;
    num_consec_best = 0;    % number of consecutive run that best solution does not improve
    num_consec_best_all = [];
    prvs_best = inf;        % previous best solution
    alpha_tt = [];          % upcrossing bound history
    beta_tt = [];           % upcrossing bound history
    Etrue_T = [];           % the transformed and normalized local minima history
    stop11 = [];
    stop2a1 = [];
    stop2b1 = [];
    fmin = inf;             % the minimal local minimum obtained so far
    fmax = -inf;            % the maximal local minimum obtained so far
    
    maxrun = 1/epsilon^2-1+((1-exp(-1))/delta-1)*(1-exp(-1))/(1-exp(-1)-delta);
    
    [cs,ct,ucs,uct] = box_con_ucon(lb,ub);
    
    BestSol.value = inf;
    BestSol.vector = unifrnd(LB,UB);
    
    if(nargin<16)
        seed_i = unifrnd(LB,UB);
    else
        seed_i = ini_seed;
    end    
    
    for i=1:10*dim
        
        if(i==1)
            pop(i).vector = seed_i;
        else
            pop(i).vector = unifrnd(LB,UB);
        end
        
        pop(i).value=objfun_i(pop(i).vector);

        maxover = maxover+1;
        
        if(~isempty(nonlcon_i))
            [c_pop,ceq_pop] = nonlcon_i(pop(i).vector);
            [c_bst,ceq_bst] = nonlcon_i(BestSol.vector);
        else
            c_pop = 0;
            ceq_pop = 0;
            c_bst = 0;
            ceq_bst = 0;
        end
        
        if isempty(ceq_pop)
            ceq_pop = 0;
        end
        if isempty(ceq_bst)
            ceq_bst = 0;
        end
        
        if pop(i).value<BestSol.value && sum(c_pop>0)==0 && sum(abs(ceq_pop)>tol_eq_i)==0 % better and feasible
            BestSol = pop(i);
            CurrentBest = pop(i);            
        elseif (pop(i).value + sum(max(0,c_pop)) + ceq_pop * ceq_pop') <...
                    (BestSol.value + sum(max(0,c_bst)) + ceq_bst * ceq_bst')       % better fittness and infeasible
            BestSol = pop(i);            
            CurrentBest = pop(i);
            CurrentBest.value = inf;            
        else                                    % worse and infeasible
            BestSol = pop(i);
            BestSol.value = inf;
            CurrentBest = pop(i);
            CurrentBest.value = inf;            
        end        
    end
    pop = [pop CurrentBest];        
    
    while(flag~=1)        
        if(num==1)
            x_ini = rand(1,dim);                    % generate the extremum seed with Gibbs sampling        
        end
        
        [xo,eo,maxover] = local_m(x_ini);   % solve the optimization problem with x_ini and obtain the design vector and function value    
        
        if(maxover==1)
            break;
        end
                
        x_i = store(x_ini,xo,eo,1);                 % store this local minimum into the observation set           
        x_nb = find_neighbor(ObSet(end));         % find the closest extremum seed to x_ini
        
        sfg = check_stop(eo,epsilon,delta);
        stopflg = [stopflg sfg];        
        
        if(stopflg(end)==1)
            flag=1;            
        end
        
        if(norm((x_i.vector-x_nb.vector).*cs)/sqrt(dim)<tol1&&abs(x_i.value-x_nb.value)<tol2)
        else
            xt = iterative_tunneling(x_i,x_nb,1);   % iterative tunneling starts from the first level        
        end          
        x_ini = Sampling();                         % generate initial points
        x_inii = [x_inii;x_ini];
        
        if(stopflg(end)==1)            
            flag=1;            
        end        
        
    end
    z1 = CurrentBest;    
    z2 = Numfeval;
end