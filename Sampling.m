function z1 = Sampling()
% this is the main function of the sampling method. z1 is the returned new
% initial point.

global dim LB UB CurrentBest BestSol maxover pop cs ct objfun_i nonlcon_i tol_eq_i
% global ucs uct

nVar=dim;           % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin= LB;         % Lower Bound of Decision Variables
VarMax= UB;         % Upper Bound of Decision Variables

beta_min=0;         % Lower Bound of Scaling Factor
beta_max=1;         % Upper Bound of Scaling Factor

pCR=0.2;            % Crossover Probability

% Sampling Main Loop
for i=1:length(pop)
    
    x = pop(i).vector;
    
    A = randperm(10*dim);
    
    A(A==i) = [];
    
    a=A(1);
    b=A(2);
    c=A(3);
    
    % Mutation
    beta=unifrnd(beta_min,beta_max,VarSize);
    y = pop(a).vector+beta.*(pop(b).vector-pop(c).vector);
    y = max(y, VarMin);
    y = min(y, VarMax);
    
    % Crossover
    z=zeros(size(x));
    j0=randi([1 numel(x)]);
    for j=1:numel(x)
        if j==j0 || rand<=pCR
            z(j)=y(j);
        else
            z(j)=x(j);
        end
    end
    
    NewSol.vector=z;
    NewSol.value=objfun_i(NewSol.vector);
    maxover = maxover+1;
    
    if(~isempty(nonlcon_i))
        [c_new,ceq_new] = nonlcon_i(NewSol.vector);
        [c_pop,ceq_pop] = nonlcon_i(pop(i).vector);
    else
        c_new = 0;
        ceq_new = 0;
        c_pop = 0;
        ceq_pop = 0;
    end
    
    if isempty(ceq_new)
        ceq_new = 0;
    end
    if isempty(ceq_pop)
        ceq_pop = 0;
    end
    
    if (NewSol.value + sum(max(0,c_new)) + ceq_new * ceq_new') <...
            (pop(i).value + sum(max(0,c_pop)) + ceq_pop * ceq_pop')        % better fitness
        pop(i)=NewSol;
        if(~isempty(nonlcon_i))
            [c_bst,ceq_bst] = nonlcon_i(BestSol.vector);
        else
            c_bst = 0;
            ceq_bst = 0;
        end
        
        if isempty(ceq_bst)
            ceq_bst = 0;
        end
        
        if (pop(i).value + sum(max(0,c_pop)) + ceq_pop * ceq_pop') <...
                (BestSol.value + sum(max(0,c_bst)) + ceq_bst * ceq_bst')
            BestSol=pop(i);
        end
        
        if(~isempty(nonlcon_i))
            [c_pop,ceq_pop] = nonlcon_i(pop(i).vector);
        else
            c_pop = 0;
            ceq_pop = 0;
        end
        
        if isempty(ceq_pop)
            ceq_pop = 0;
        end
        if pop(i).value<CurrentBest.value && sum(c_pop>0)==0 && sum(abs(ceq_pop)>tol_eq_i)==0
            CurrentBest=pop(i);
        end
    end
end

% Update Best Cost
z1 = BestSol.vector.*cs+ct;
end