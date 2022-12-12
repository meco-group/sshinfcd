classdef postprocessor 
% POSTPROCESSOR Postprocessing steps to the controller matrices originating from the standard problem

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
    
    properties
        % prob - The hinfcd problem
        % The extended H-infinity controller synthesis problem.
        prob
        
        % K - The controller resulting from the standard problem
        K
        
        % perf - Achieved performance
        perf
    end
    
    methods
        %% Constructor
        function obj = postprocessor(problem, K, perf)
        % POSTPROCESSOR Constructor
            
            % construct the object            
            assert(isa(problem,'sshinfcd.problem'), 'Invalid problem.'); 
            obj.prob = problem; 
            obj.K = K; 
            obj.perf = perf;
            
        end
        
        %% Recombination the controller
        function obj = recombine(obj)
        % RECOMBINE Adds the impulsive and unstable modes of the weighting filters to the controller
            
            % controller parameters from the standard problem
            Ac = obj.K.A; 
            Bc1 = obj.K.B(:,1:obj.prob.no()); Bc2 = obj.K.B(:,(obj.prob.no()+1):end);
            Cc1 = obj.K.C(1:obj.prob.ni(),:); Cc2 = obj.K.C((obj.prob.ni()+1):end,:); Cc = [Cc1 ; Cc2];
            Dc11 = obj.K.D(1:obj.prob.ni(),1:obj.prob.no()); Dc12 = obj.K.D(1:obj.prob.ni(),(obj.prob.no()+1):end);
            Dc21 = obj.K.D((obj.prob.ni()+1):end,1:obj.prob.no()); Dc22 = obj.K.D((obj.prob.ni()+1):end,(obj.prob.no()+1):end);
            Dc1 = [Dc11 ; Dc21]; Dc2 = [Dc12 ; Dc22];
            
            % separation variables
            [GAMMAo,~,Zotilde] = ofseparator(obj.prob); 
            [GAMMAi,~,Zi] = ifseparator(obj.prob);
            
            % reconstruct controller parameters
            A = [obj.prob.Ao()-Zotilde*Dc1      -Zotilde*Cc     -(GAMMAo-Zotilde*Dc2)*Zi;
                 Bc1                            Ac              -Bc2*Zi                 ;
                 Dc11                           Cc1             obj.prob.Ai()-Dc12*Zi   ];
            B = [GAMMAo-Zotilde*Dc2;
                 Bc2               ;
                 Dc12              ];
            C = [Dc21 Cc2 GAMMAi-Dc22*Zi];
            D = Dc22;
            
            % return 
            obj.K = ss(A,B,C,D);
        end
        
        %% Compensation for the feedthrough component of the generalized plant
        function obj = compensateDyu(obj)
        % COMPENSATEDYU Compensates for direct feedthrough in the generalized plant, if present
            
            % get synthesized controller
            [A0,B0,C0,D0] = dssdata(obj.K); 
            
            % get feedthrough matrix
            Dyu = obj.prob.Dyu(); 
            
            % compensate for the algebraic loop
            Mcyu = eye(size(obj.K,1))+D0*Dyu;
            C = Mcyu\C0;
            D = Mcyu\D0;
            A = A0-B0*Dyu*C;
            B = B0*(eye(size(obj.K,2))-Dyu*D); 
            
            % set updated controller data
            obj.K = set(obj.K,'A',A,'B',B,'C',C,'D',D);
            
        end
           
    end
    
end