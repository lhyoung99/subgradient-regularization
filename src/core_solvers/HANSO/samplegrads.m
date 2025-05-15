function [X,G,oracle] = samplegrads(x,eps,ngrad,func_value,subgrad)
%
% to be called only by gradsamp
% sample gradients at points generated uniformly in the 
% 2-norm ball of radius eps around x
% (in the original 2005 implementation we used the inf-norm ball)

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

nvar = length(x);
X = zeros(nvar,ngrad);
G = zeros(nvar,ngrad);
oracle = 0;
for j = 1:ngrad 
   % uniform distribution over inf-norm ball
   % xpert = x + 2*eps*(rand(nvar,1) - 0.5); 
   %
   % uniform distribution over 2-norm ball: from Frank Curtis
   u = randn(nvar,1); % note: randn
   u = eps*(rand^(1/nvar))*u/norm(u); % note: rand
   xpert = x + u;
   % fprintf('%.4f %.4f\n', xpert);
   fpert = func_value(xpert);
   gpert = subgrad(xpert); 
   oracle = oracle + 1;
   count = 0;
   while isnaninf(fpert) || isnaninf(gpert)  % in particular, disallow infinite function values
       count = count + 1;
       if count > 100 % should never happen, but just in case
           error('gradsamp: too many contractions needed to find finite f and grad values')
       end
       xpert = (x + xpert)/2;     % contract back until feasible
       % fpert = func_value(xpert);
       gpert = subgrad(xpert); 
       oracle = oracle + 1;
   end % discard function values
   X(:,j) = xpert;
   G(:,j) = gpert;   
end