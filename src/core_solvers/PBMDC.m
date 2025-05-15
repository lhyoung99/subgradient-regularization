function results = PBMDC(...
    x0, ... 
    func_val1, ...
    func_val2, ...
    subgrad1, ...
    subgrad2, ...
    flags, ... 
    method_params ...
)
%--------------------------------------------------------------------------
%         Proximal bundle method for DC programming 
%
%               By Welington de Oliveira
%
% PBMDC computes a critical point for DC programs of the form
%
%            min  f1(x) - f2(x)    s.t. x    \in X
% 
% where f1,f2: R^n \to R  are convex functions and X is either a 
% polyhedron or the whole space R^n
%
% The problem to be solved must be defined in function PbmData.m
%
%--------------------------------------------------------------------------
% This is an implemetation of Algorithm 1 and 2 given in the manuscript
%
% [1] "Proximal bundle methods for nonsmooth DC programming", June 2018
%  by Welington de Oliveira
%
% The manuscript is available at: http://www.oliveira.mat.br/publications
%
% If f2 is the pointwise maximum of finitely many functions and the
% provided oracles returns all the gradients of the active functions, then
% this function is ensured to compute a d-stationary point of the problem.
% See [1] for more details. 
% See the function pbm11ctive.m for an example of
% an oracle providing gradients of the active functions
% 
% The algorithm uses the simpler descent test (12) given in [1] and employs
% Matlab quadprog for solving the QP master problem.
% In [1] Gurobi was used instead of quadprog. For that, set 
% method_params.matlab = false in the function PbmData.m
% 
%--------------------------------------------------------------------------
% If you use this matlab function, please cite [1] in your work.
% If you find any bug, please contact me www.oliveira.mat.br
%--------------------------------------------------------------------------
%
%=========================================================================%
% [Modified in 2025]
% Returns:
%     A struct containing hisory results
%=========================================================================%

% Initialization
method_name = 'PBMDC';
% tic;
ns = 0;ss=0;results.ind=-1;
%
step_tol = method_params.step_tol;
obj_tol = method_params.obj_tol;
max_oracle_call = method_params.max_oracle_call;
max_time_seconds = method_params.max_time_seconds;
print_frequency = method_params.print_frequency;

opt_val = [];
if isfield(flags, 'opt_val')
    opt_val = flags.opt_val;
end

opt_soln = [];
if isfield(flags, 'opt_soln')
    opt_soln = flags.opt_soln;
end

% Start timer
start_time = tic;

% First stability center
xc = x0;
f1 = func_val1(xc); f2 = func_val2(xc);
fxc = f1 - f2;
obj_best = fxc;

g1 = subgrad1(xc); g2 = subgrad2(xc);
num_oracle = 1;

% Record initial objective and distance
if ~isempty(opt_soln)
    distance_current = norm(xc - opt_soln);
else
    distance_current = [];
end
objective_values_list = fxc; 
history = xc(:);
distance_list = distance_or_empty(distance_current);
num_oracle_list = num_oracle;
elapsed_time_list = toc(start_time);

%--------------------------------------------------------------------------
% Cutting plane of f1 : f1(x) + g1'(y-x) \geq r ==> g1'y - r \geq g1'x - f1(x)
modeloPC(1).G   = g1';
modeloPC(1).rhs = g1'*xc - f1;
% Cutting plane of f2 : f2(x) + g2'(y-x) \geq r ==> g2'y - r \geq g2'x - f2(x))
nactive         = size(g2,2);
modeloPC(2).G   = g2';
modeloPC(2).rhs = g2'*xc - f2*ones(nactive,1);            
%--------------------------------------------------------------------------
% Loop
k = 0;nbun=1;nbun2=1;in=0;ds=0;
while true
    k = k + 1;
    % Computing a new point x_{k+1}
    if method_params.QPdual
       [x_new,mu,fcheck] = SubProbProximalDualDC(xc,modeloPC,method_params);
    else
       [x_new,mu,fcheck] = SubProbProximalDC(xc,modeloPC,method_params);
    end
    
    if ~isempty(opt_soln)
        distance_current = norm(x_new - opt_soln);
    end

    % Oracle call        
    % [f1,g1,f2,g2,ind] = feval(oracle,pars);
    % if ind < 0; fprintf(1,'%s',results.msg); return;end
    f1 = func_val1(x_new); f2 = func_val2(x_new);
    g1 = subgrad1(x_new); g2 = subgrad2(x_new);
    obj_best = min([obj_best, f1 - f2]);
    num_oracle = num_oracle + 1;
    
    % Record objective and distance
    objective_values_list(end+1) = obj_best; 
    history = [history, x_new(:)];
    distance_list(end+1) = distance_or_empty(distance_current);
    num_oracle_list(end+1) = num_oracle;
    elapsed_time_list(end+1) = toc(start_time);

    %----------------------------------------------------------------------       
    % Stopping test    
    normdif = norm(x_new-xc);
    %if normdif <=step_tol
    if should_terminate(...
            elapsed_time_list(end), ...
            max_time_seconds, ...
            num_oracle, ...
            max_oracle_call, ...
            f1 - f2, ...
            opt_val, ...
            obj_tol, ...
            normdif, ...
            step_tol)
        % Print final info
        fprintf(1,'k: %3d, |xkk -xck|: %4.1d  fup: %4.1d \n',k,normdif,obj_best);

        results = compile_statistics(method_name, x_new, num_oracle_list, ...
                elapsed_time_list, objective_values_list, distance_list, history, flags);
        return;
    end
    %----------------------------------------------------------------------       

    % Descent test    
    aux = fxc - 0.5*(method_params.kappa/method_params.tmax)*normdif^2;    
    ps=0;
    if f1-f2 <= aux  
        % Serious step
        ds      = ds+1;         
        xc = x_new;
        fxc     = f1 - f2;
        ss      = ss+1;

        

        if ss>3, method_params.t = min(method_params.tmax,method_params.t*(1.25));ss=0;end
        ps = 1; ns = 0; in = 0;
        % Updating the cutting-plane model of f2, if a DC model is used
        % If a convex model is employed, then method_params.nBun2 = 1 and only one
        % linearization of f2 is considered
        nactive = size(g2,2);
        if nactive>1
            modeloPC(2).G   = g2';
            modeloPC(2).rhs = g2'*xc - f2*ones(nactive,1);            
        elseif(nbun2<method_params.nBun2)
            modeloPC(2).G   = [modeloPC(2).G; g2'];
            modeloPC(2).rhs = [modeloPC(2).rhs; g2'*xc - f2];       
            nbun2           = nbun2 + 1;
        else
            modeloPC(2).G   = [modeloPC(2).G(2:end,:); g2'];
            modeloPC(2).rhs = [modeloPC(2).rhs(2:end); g2'*xc - f2];       
            nbun2           = method_params.nBun2;
        end
    else
        % Null step
        ns = ns+1;
        ss = 0;
         if ns>2
             method_params.t=max(method_params.tmin,method_params.t/(in+2));
             in = in+1; ns=0;
         end
    end
   %----------------------------------------------------------------------
   % Display progress on the screen
   if mod(k, print_frequency) == 0
       fprintf(1,'k: %3d, |xkk -xck|: %4.1d  fup: %8.4d, t %4.2d ,  nb %3d,  nb2 %3d,  ps %3d  \n',...
                  k,       normdif,   fxc, method_params.t, nbun, nbun2, ps );
   end
   %----------------------------------------------------------------------

   % Managing the bundle of information of f1
   J  = (abs(mu) > 1e-7);
   nAc=sum(J);nNAc = nbun-nAc;
   if (nbun>10*method_params.nBun);ps=1;end
   if (ps==0)||(nbun<method_params.nBun)
      modeloPC(1).G   = [modeloPC(1).G; g1'];
      modeloPC(1).rhs = [modeloPC(1).rhs; g1'*x_new - f1];       
   elseif (nAc>=method_params.nBun)
       i = nAc-method_params.nBun;
       aux = [1:nbun];J=aux(J);
       if i>3
           J = J(i+3:end);
       else 
           J=J(3:end);
       end
       d = modeloPC(1).G'*mu;
       modeloPC(1).G   = [modeloPC(1).G(J,:);d'; g1'];
       modeloPC(1).rhs = [modeloPC(1).rhs(J);d'*x_new - fcheck; g1'*x_new - f1];
   else  
       i = method_params.nBun-nAc -1;
       aux = 1:nbun;I=aux(~J);I=I(nNAc-i+1:nNAc);
       modeloPC(1).G   = [modeloPC(1).G(I,:);modeloPC(1).G(J,:);g1'];
       modeloPC(1).rhs = [modeloPC(1).rhs(I);modeloPC(1).rhs(J); g1'*x_new - f1];
   end
   nbun = length(modeloPC(1).rhs);
end
%==========================================================================
% End of processing
% if k > method_params.max_iter 
%     results.msg = 'Maximum number of iterations';
% else
%     results.ind=1; 
%     results.msg = 'Criticality!'; 
% end
% results.k   = k;        % number of iterations
% results.sol = xc;  % Best computed point
% results.val = fxc;      % Best value of the function
% results.cpu = toc;      % CPU time
% results.ds  = ds;       % Number of serious steps
% return

end


% #########################################################################
%                QP Master problem
% 
% #########################################################################
function [x_new,Dualvar,fcheck] = SubProbProximalDualDC(xc, ...
                                modeloPC, ...
                                method_params)
%
% This function solves the master program defining the next trial point
% by its dual. This is only possible if the DC problem has no constraint
%
% If a DC model is employed (ncut2>1), then this function computes the trial
% point by globally solving the master program
%
% 
ncut1 = length(modeloPC(1).rhs);
ncut2 = length(modeloPC(2).rhs);
Aeq   = ones(1,ncut1);
beq   = 1;
lb    = zeros(ncut1,1);
model = inf;
for i=1:ncut2    
   K     = (repmat(modeloPC(2).G(i,:),ncut1,1) - modeloPC(1).G  )';
   Q     = method_params.t*(K'*K);
   condK = 1./condest(Q);
   tr    = 0; 
   sigma = 1.d-12;
   while condK < 1e-15 && tr<= 10
     tr    = tr + 1;
     Q     = Q + sigma*eye(ncut1);
     sigma = 10*sigma;
     condK = 1./condest(Q);
   end
   q   = modeloPC(1).rhs + K'*xc;
   %-------------------------------------------------------------------------
   if method_params.matlab
     % Solving the QP by quadprog: NOT RELIABLE!
     % opts = optimset('Display','off');
     opts.Display = 'off'; 
     % opts.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
     % mu0 = zeros(ncut1,1);
     % [~,idx] = min(q);
     % mu0(idx) = 1;
     [mu,~,ind]  = quadprog(Q,q,[],[],Aeq,beq,lb,[],[],opts);
   else
     % Solving the QP by Gurobi
     prob.obj   = q;
     prob.Q     = 0.5*sparse(Q);
     prob.A     = sparse(Aeq);
     prob.rhs   = beq;
     prob.sense = '=';
     prob.lb    = lb;
     % par.resultputFlag =0; % 1 to print resultput
     par.OutputFlag = 0;
     res = gurobi(prob,par);
     if strcmp(res.status, 'OPTIMAL')|| strcmp(res.status, 'SUBOPTIMAL')
        mu   = res.x;
        ind = 1;
     else
       ind=-1;
     end
   end
   if ind>=0
     x       = xc + method_params.t*K*mu;
     fcheck1 = max(modeloPC(1).G*x - modeloPC(1).rhs);
     fcheck2 = max(modeloPC(2).G*x - modeloPC(2).rhs);
     %
     if fcheck1-fcheck2< model
        model  = fcheck1-fcheck2;
        fcheck = fcheck1;
        x_new = x;
        Dualvar= mu;
     end
   end
end
return
end

% ########################################################################
function [x_new,Dualvar,fcheck] = SubProbProximalDC(xc, ...
                                modeloPC, ...
                                method_params)
%
% This function solves the master program defining the next trial point
% by its primal.
%
% If a DC model is employed (ncut2>1), then this function computes the trial
% point by globally solving the master program
%
n    = length(xc);                     
ncut1 = length(modeloPC(1).rhs);   
ncut2 = length(modeloPC(2).rhs);   
model = inf;
%-------------------------------------------------------------------------
% Add the equality constraints (if any) of the problem
% Aeq = method_params.X.Aeq;
% beq = method_params.X.beq;
% neq = length(beq);
% if ~isempty(Aeq)
%     m   = size(Aeq,1);
%     Aeq = [Aeq, zeros(m,1)];
% end
%-------------------------------------------------------------------------
% Add the inequality constraints (if any) of the problem
% A = method_params.X.A;
% b = method_params.X.b;
% nineq = length(b);
% if ~isempty(A)
%     m   = size(A,1);
%     A   = [A,        zeros(m,1);
%            modeloPC(1).G, -ones(ncut1,1)];
%     b   = [b;modeloPC(1).rhs];
% else
%     A   = [modeloPC(1).G, -ones(ncut1,1)];
%     b   = [modeloPC(1).rhs];
% end
%-------------------------------------------------------------------------
% Add bounds (if any) of the problem
% lb = method_params.X.lb;
% if ~isempty(lb)
%     lb = [lb;method_params.flow];
% else
%     lb      = -inf(n+1,1);
%     lb(n+1) = method_params.flow;
% end
% ub = method_params.X.ub;
% if ~isempty(ub)
%     ub = [ub;inf];
% else
%     ub      = inf(n+1,1);
% end
%-------------------------------------------------------------------------
% Define the objective function
aux       = 1/method_params.t;
Q         = diag(ones(n+1,1)*aux);
Q(n+1,n+1)= 0;
% m         = length(b);
prob.Q    = 0.5*sparse(Q);
% prob.A    = sparse([A;Aeq]);
% prob.rhs  = [b;beq];
% m1        = length(prob.rhs);
% sense=[];for i=1:m;sense=strcat(sense,'<');end
% for j=m+1:m1; sense=strcat(sense,'=');end
% prob.sense = sense;
% prob.lb    = lb;
% prob.ub    = ub;
prob.lb = -inf(n+1,1);
Dualvar    = ones(ncut1,1)/ncut1; fcheck = -inf;
% par.resultputFlag =0; % 1 to print resultput
par.OutputFlag = 0;
for i=1:ncut2
    q =  [-modeloPC(2).G(i,:)'-aux*xc;1];
    if method_params.matlab
        % Solving the QP by quadprog
        % opts = optimset('Display','off');
        opts.Display = 'off';
        % opts.Algorithm = 'active-set';  % 'interior-point-convex' (default)/'trust-region-reflective'/'active-set'
        [x,~,ind,~,lambdap]  = quadprog(Q,q,A,b,Aeq,beq,lb,ub,[],opts);
        mu                         = lambdap.ineqlin;
        mu                         = mu(nineq+neq+1:end);
    else
       % Solving the QP by Gurobi
        prob.obj = q;
        res = gurobi(prob,par);
        if strcmp(res.status, 'OPTIMAL')|| strcmp(res.status, 'SUBOPTIMAL')
            x=res.x;
            ind = 1;
            mu  = -res.pi;
            mu  = mu(nineq+neq+1:end);
        else
            ind=-1;
        end
    end
   if ind>=0
     fcheck1  = x(end);
     x       = x(1:end-1);
     fcheck2 = max(modeloPC(2).G*x - modeloPC(2).rhs);
     if fcheck1-fcheck2< model
        model  = fcheck1-fcheck2;
        fcheck = fcheck1;
        x_new = x;
        Dualvar= mu;
     end
   end
end
if fcheck==-inf
    error('Numerical issues when solving the master program');
end
return
end

% [Added in 2025]
%% ===================== Helper functions =====================
function val = distance_or_empty(distance_current)
    if isempty(distance_current)
        val = NaN; 
    else
        val = distance_current;
    end
end

function stop_flag = should_terminate(...
                    elapsed_time, ...
                    max_time_seconds, ...
                    num_oracle, ...
                    max_oracle_call, ...
                    obj_current, ...
                    opt_val, ...
                    obj_tol, ...
                    normdif, ...
                    step_tol)
    stop_flag = false;
    if (elapsed_time > max_time_seconds) || (num_oracle >= max_oracle_call)
        disp('Maximum oracle calls or time limit reached.');
        stop_flag = true;
        return;
    end
    
    if ~isempty(opt_val) && ~isempty(obj_tol) && ((obj_current - opt_val) <= obj_tol)
        disp('Objective tolerance reached.');
        stop_flag = true;
        return;
    end
    
    if ~isempty(step_tol) && (normdif <= step_tol)
        disp('Tolerance of $\|x_{k+1} - x_k\|$ reached.');
        stop_flag = true;
        return;
    end
end

function results_struct = compile_statistics(...
    method_name, x_k, num_oracle_list, elapsed_time_list, ...
    objective_values_list, distance_list, history, flags ...
)
    results_struct.name_of_method = method_name;
    results_struct.x = x_k;
    % results_struct.num_oracle_list = 1:num_oracle;
    results_struct.num_oracle_list = num_oracle_list;
    results_struct.elapsed_time_list = elapsed_time_list;
    results_struct.objective_values_list = objective_values_list;
    results_struct.distance_list = distance_list;
    results_struct.history = history;
    results_struct.flags = flags;
end