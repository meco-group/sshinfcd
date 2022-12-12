function s = options(obj)
% OPTIONS Returns the standard settings and allows the user to easily change his/her preferred options in the structure
    
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

    % synthesis
    s.synthesis.gammasolver = 'lmilab';     % 'lmilab', 'cvx', 'yalmip'
    s.synthesis.zerotol = 0;                % any real number >= 0.
    s.synthesis.lyapunovshape = 3;          % 1 (-gamma,-gamma), 2 (-gamma^2,-1), 3 (-1,-gamma^2) or 4 (try all three)
    
    % LMILAB interface
    s.lmilab = [1e-6 500 0 0 0];           % LMILAB options, see the mincx() and feasp() documentation
    
    % YALMIP interface
    try
        s.yalmip = sdpsettings();           % YALMIP options, see the YALMIP documentation
    catch 
        s.yalmip = struct(); 
    end
    
    % CVX interface
    s.cvx.solver = '';                      % solver for CVX, see the CVX documentation 
    s.cvx.solver_settings = {''};           % solver settings for CVX, see the CVX documentation
    s.cvx.precision = 'default';            % precision settings for CVX, see the CVX documentation
    
end