function [xnew,fnew,gnew,evals,ls_termcode,fline] = linesch_gradsamp(x,f,d,func_value,subgrad,options)
%
% to be called only by gradsamp
% backtracking line search for Gradient Sampling with nonstandard Armijo condition
% see Algorithm GS in  
%    J.V. Burke, F.E. Curtis, A.S. Lewis, M.L. Overton and L.E.A. Simões, 
%    Gradient Sampling Methods for Nonsmooth Optimization,
%    In: Numerical Nonsmooth Optimization, edited by A. Bagirov et al, 
%    Springer (2020), pp. 201-225, https://arxiv.org/abs/1804.11003.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  HANSO 3.0 Copyright (C) 2021  Michael Overton
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xnorm = norm(x); % for rounding error test
dnorm = norm(d); 
beta = options.beta; % for the non-standard Armijo condition
betadnormsq = beta*dnorm^2; % for the non-standard Armijo condition 
gamma = options.gamma; % contraction factor
delta = options.delta; % rounding error termination test
t_lowbound = delta*xnorm/dnorm;
prtlevel = options.prtlevel;
evals = 0;
t = 1;
done = false;
while ~done
    xnew = x + t*d; % sign of d was changed in calling routine
    fnew = func_value(xnew);
    gnew = subgrad(xnew);
    evals = evals + 1;
    fline(evals) = fnew;
    if prtlevel > 2
        fprintf('line search: t=%g, evals=%d, fnew=%22.16e\n',t,evals,fnew)
    end
    if fnew < f - t*betadnormsq % non-standard Armijo condition
        ls_termcode = 0;
        done = true;
    elseif t < t_lowbound % if delta is machine epsilon, no sense continuing
        ls_termcode = 1;
        done = true;
    else
        t = gamma*t;
    end
end