function X = strucsylvester(A,B,C,nx)
% STRUCSYLVESTER Solves a structured Sylvester equation
%
% [X] = STRUCSYLVESTER(A,B,C,nx) returns the solution to the 
% generalized Sylvester equation with specific structure
%   A*[X1;X2] + [X1;0]*B = C
%
% with size(X1,1) = nx. 
%
% See also HINFCD.UTIL.GENSYLVESTER.

% This file is part of sshinfcd.
% Copyright (c) 2022, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% sshinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% sshinfcd is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
% License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with sshinfcd. If not, see <https://www.gnu.org/licenses/>.

    % get dimensions
    m = size(B,1); n = size(A,1); r = size(A,2);
    
    % solve problem
    UL = kron(eye(m),A);
    UR = kron(B', [[eye(nx) ; zeros(n-nx,nx)] zeros(n,r-nx)]);
    sol = (UL+UR)\C(:);
    
    % check relative residual
    res = (UL+UR)*sol-C(:);
    assert(isempty(res) || max(abs(res)) <= sqrt(eps), 'STRUCSYLVESTER: Absolute residual is bigger than sqrt(eps).'); 
    
    % return solution
    X = reshape(sol,[r,m]);
    
end