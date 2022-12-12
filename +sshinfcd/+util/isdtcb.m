function d = isdtcb(G)
% ISDTCB Checks detectability of an LTI system

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

    [A,~,~,~,k] = obsvf(G.A,G.B,G.C);
    k = length(A)-sum(k);
    d = all(real(eig(A(1:k,1:k)))<0); 
end
