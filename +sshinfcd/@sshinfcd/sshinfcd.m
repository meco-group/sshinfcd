classdef sshinfcd
% SSHINFCD Tool for state-space H-infinity controller synthesis through LMIs

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

    properties (Access=public)
        % problem - The control problem
        % The control problem contains the generalized plant, input and
        % output weighting filters and the specifications.
        % See also SSHINFCD.PROBLEM.
        problem
        
        % opts - Options
        % The options allow the user to specify his own preferences for
        % selecting solvers, tolerances, etc. 
        % See also SSHINFCD.OPTIONS.
        opts
        
        % preprocessor - The problem preprocessor
        % The problem preprocessor transforms the problem into a standard
        % formulation that allows reformulation into an SDP. 
        % See also SSHINFCD.PREPROCESSOR.
        preprocessor
              
        % postprocessor - The problem postprocessor
        % The problem postprocessor combines the solution of the SDP and 
        % and the information of the preprocessing phase into the
        % controller that is optimized for. 
        % See also SSHINFCD.POSTPROCESSOR.
        postprocessor
        
        % watchdog - Monitor of the solution progress
        % The watchdog keeps track of the evolution of the different stages 
        % of the solution, save this information in a logbook, reports to
        % the user in the command window, and is responsible for error 
        % handling.  
        % See also SSHINFCD.WATCHDOG.
        watchdog = sshinfcd.watchdog.watchdog(); 
    end
    
    methods
        function obj = sshinfcd()
        end
        
        function obj = setproblem(obj,P,Wi,Wo,specs,region)
        % SETPROBLEM Initializes the control problem
            if nargin < 6 % no pole region constraints
                region = struct('L',[],'M',[]); 
            end
            obj.problem = sshinfcd.problem(obj.watchdog,P,Wi,Wo,specs,region);
        end
        
        obj = setoptions(obj,s);
        obj = options(obj);
        [obj,gamma] = solve(obj,options);  
    end
    
end
