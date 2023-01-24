classdef preprocessor 
% PREPROCESSOR Processor formulating a problem into an SDP

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
        % prob - The sshinfcd problem
        % The extended H-infinity controller synthesis problem.
        prob
        
        % opts - The options for controller synthesis
        opts
    end
    
    methods
        %% Constructor
        function obj = preprocessor(problem, options)
        % PREPROCESSOR Constructor
            
            % construct the object            
            assert(isa(problem,'sshinfcd.problem'), 'Invalid problem.'); 
            obj.prob = problem; 
            
            assert(isstruct(options), 'options should be a structure.');
            obj.opts = options;
            
        end
        
        %% Solution procedure
        function [K,gamma] = solve(obj)
            % SOLVE Executes the solution routines
            
                iscstr = [obj.prob.specs.weight]==0;
                isobj = ~iscstr;
                opts_ = obj.opts; 
                
                % Step 1: Solve a first time to estimate the scaling of the problem
                    if obj.opts.synthesis.lyapunovshape == 4
                        opts_.synthesis.lyapunovshape = 1; obj1 = Inf;
                        try
                            [~,gamma{1}] = sshinfcd.standard.nlcov.nlcov(obj.prob, opts_);
                            obj1 = obj.prob.specs(isobj).weight*gamma{1};
                        catch
                            % do nothing
                        end
                        opts_.synthesis.lyapunovshape = 2; obj2 = Inf;
                        try
                            [~,gamma{2}] = sshinfcd.standard.nlcov.nlcov(obj.prob, opts_);
                            obj2 = obj.prob.specs(isobj).weight*gamma{2};
                        catch
                            % do nothing
                        end
                        opts_.synthesis.lyapunovshape = 3; obj3 = Inf; 
                        try
                            [~,gamma{3}] = sshinfcd.standard.nlcov.nlcov(obj.prob, opts_);
                            obj3 = obj.prob.specs(isobj).weight*gamma{3};
                        catch
                            % do nothing
                        end
                        [~,best] = min([obj1 obj2 obj3]);
                        gamma = gamma{best}; 
                        opts_.synthesis.lyapunovshape = best; 
                    else
                        [~,gamma] = sshinfcd.standard.nlcov.nlcov(obj.prob, opts_); 
                    end
                    if opts_.synthesis.lyapunovshape == 2 || opts_.synthesis.lyapunovshape == 3
                        gamma = sqrt(gamma);
                    end
                    
                % Step 2: Rescale the problem 
                PRE = eye(obj.prob.nw());
                POST = eye(obj.prob.nz());
                switch opts_.synthesis.lyapunovshape 
                    case 1
                        for i=1:length(obj.prob.specs)
                            j = 1; 
                            if isobj(i)
                                PRE(obj.prob.specs(i).in,obj.prob.specs(i).in) = 1/sqrt(gamma(j))*eye(length(obj.prob.specs(i).in));
                                POST(obj.prob.specs(i).out,obj.prob.specs(i).out) = 1/sqrt(gamma(j))*eye(length(obj.prob.specs(i).in));
                                obj.prob.specs(i).weight = obj.prob.specs(i).weight*gamma(j); 
                                j = j+1;
                            end
                        end
                    case 2
                        for i=1:length(obj.prob.specs)
                            j = 1; 
                            if isobj(i)
                                PRE(obj.prob.specs(i).in,obj.prob.specs(i).in) = 1/gamma(j)*eye(length(obj.prob.specs(i).in));
                                obj.prob.specs(i).weight = obj.prob.specs(i).weight*gamma(j)^2; 
                                j = j+1;
                            end
                        end
                    case 3
                        for i=1:length(obj.prob.specs)
                            j = 1; 
                            if isobj(i)
                                POST(obj.prob.specs(i).out,obj.prob.specs(i).out) = 1/gamma(j)*eye(length(obj.prob.specs(i).out));
                                obj.prob.specs(i).weight = obj.prob.specs(i).weight*gamma(j)^2; 
                                j = j+1;
                            end
                        end
                end
                obj.prob.Wi = obj.prob.Wi*PRE;
                obj.prob.Wo = POST*obj.prob.Wo; 
                
                % Step 3: Solve again with (hopefully) less conservatism
                % Ideally, one would have to repeat this iteratively until
                % all gamma values that are optimized converge to 1. 
                
                [K,delta_gamma] = sshinfcd.standard.nlcov.nlcov(obj.prob, opts_);
                if opts_.synthesis.lyapunovshape == 2 || opts_.synthesis.lyapunovshape == 3
                    delta_gamma = sqrt(delta_gamma);
                end
                gamma = delta_gamma.*gamma; 
        end
        
    end
    
end