function loc = getloc(x, X, G, d, w)
% to be called only by hanso 
% 
% Get local optimality certificate from the set of gradients computed by
% BFGS or gradient sampling. Verify that x is the first column of X,
% and that d = G*w. These should always be the case, but costs nothing to
% check. Set loc.evaldist to the max distance from x to the other columns
% of X and set loc.dnorm to ||d||.
%

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

% double check that d = G*w within rounding error and that w>=0
dnorm = norm(d);
if any(w < 0) || norm(d - G*w) > 0
    error('getloc: w is not >= 0 or d is not equal to G*w: email mo1@nyu.edu')
end
loc.dnorm = dnorm;
% compute distances from x to other columns of X
for j = 1:size(X,2)
    dist(j) = norm(x - X(:,j)); % 
end
% double check that x equals the first column of X 
if dist(1) ~= 0
    error('getloc: x does not equal the first column of X: email mo1@nyu.edu')
end
loc.evaldist = max(dist); 