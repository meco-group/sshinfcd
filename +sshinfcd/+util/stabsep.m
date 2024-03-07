function [GS,GNS] = stabsep(G,region)
% STABSEP Separates the D-stable and the D-unstable dynamics (D described by 'region')
% Unstable modes (positive or zero real parts) are always classified as
% D-unstable, as these would be incompatible with H-infinity performance. 
% See also STABSEP. 

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

    % get the intersection of the LMI regions
    M = 1;
    L = 0;
    for i=1:length(region)
        M = blkdiag(M,region(i).M);
        L = blkdiag(L,region(i).L);
    end

    % go to the modal realization and separate modes
    G = canon(G,'modal');
    i = 1;
    is = [];
    ius = [];
    while i <= length(G.A)
        if i < length(G.A) && G.A(i+1,i) ~= 0 
            z = eig(G.A(i:i+1,i:i+1));
            if all(eig(L+M*z(1)+M'*conj(z(1)))<0) && all(eig(L+M*z(2)+M'*conj(z(2)))<0) % D-stable
                is = [is i i+1];
            else % D-unstable
                ius = [ius i i+1];
            end
            i = i+2;
        else
            if all(eig(L+M*G.A(i,i)+M'*G.A(i,i))<0) % D-stable
                is = [is i];
            else % D-unstable
                ius = [ius i];
            end
            i = i+1;
        end
    end

    if isempty(G.E); G.E = eye(size(G.A)); end
    GS = ss(G.A(is,is),G.B(is,:),G.C(:,is),G.D,G.Ts);
    GNS = ss(G.A(ius,ius),G.B(ius,:),G.C(:,ius),zeros(size(G.D)),G.Ts);
    
end