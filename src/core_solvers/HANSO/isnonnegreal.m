function ipr = isnonnegreal(x)
% return true if x is a nonnegative real scalar, false otherwise 

if ~isscalar(x) 
    ipr = 0;
else % following is OK since x is scalar
    ipr = isreal(x) & x >= 0;
end