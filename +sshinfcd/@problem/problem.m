classdef problem 
% PROBLEM Extended H-infinity controller synthesis problem

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
        % P - Generalized plant
        % A state-space realization of the (unweighted) generalized plant. 
        P

        % Wi - Input weighting filter
        % A state-space realization of the input weighting filter. 
        Wi

        % Wo - Output weighting filter
        % A state-space realization of the output weighting filter. 
        Wo

        % specs - Desired performance specifications
        % Structure array indicating the performance channels and their relative weights.
        specs
        
        % region - Desired closed-loop D-stability specificions
        % Structure array indicating the closed-loop pole region constraints.
        region

        % watchdog - A watchdog for monitoring
        % The watchdog keeps track of all processes that happen on objects of this class. 
        watchdog
        
    end
    
    methods
        %% Constructor, input checking, standard form
        
        function obj = problem(watchdog, P, Wi, Wo, specs, region)
        % PROBLEM Constructor
        
            % construct the object
            assert(isa(watchdog,'sshinfcd.watchdog.watchdog'), 'Invalid watchdog.');
            obj.watchdog = watchdog; 
            obj.P = P; 
            obj.Wi = Wi;
            obj.Wo = Wo; 
            obj.specs = specs;
            obj.region = region; 
            
            % check the validity of the input
            obj.checkInput(); 
            
            % preprocess the object
            obj = obj.standardfilters(); 
            obj = obj.balreal(); 
            obj.checkassumptions();
            obj = obj.sortchannels(); 
        end

        function checkInput(obj)
        % CHECKINPUT Checks the validity of the input arguments
        
            % check member types
            obj.watchdog.assertType(obj.watchdog, 'sshinfcd.watchdog.watchdog', 'Invalid watchdog type.');
            obj.watchdog.assertType(obj.P, 'numlti', 'P is not a standard MATLAB representation of an LTI system.'); 
            obj.watchdog.assertType(obj.Wi, 'numlti', 'Wi is not a standard MATLAB representation of an LTI system.'); 
            obj.watchdog.assertType(obj.Wo, 'numlti', 'Wo is not a standard MATLAB representation of an LTI system.');
            obj.watchdog.assertType(obj.specs, 'struct', 'specs is not a structure.'); 
            obj.watchdog.assertType(obj.region, 'struct', 'region is not a structure.'); 
            
            % check sample times
            obj.watchdog.assertCond(obj.P.Ts==obj.Wi.Ts && obj.P.Ts==obj.Wo.Ts, 'Sampling times of Wi, Wo and P should be equal.'); 

            % assert compatible sizes
            obj.watchdog.assertCond(size(obj.P,1)>size(obj.Wo,2), 'At least one measured output is required.');
            obj.watchdog.assertCond(size(obj.P,2)>size(obj.Wi,1), 'At least one control input is required.');
            
            % check specifications
            obj.watchdog.assertField(obj.specs, 'in', 'specs requires a field ''in''.');
            obj.watchdog.assertField(obj.specs, 'out', 'specs requires a field ''out''.'); 
            obj.watchdog.assertField(obj.specs, 'weight', 'specs requires a field ''weight''.'); 
            cellfun(@(x) obj.watchdog.assertNonEmpty(x, 'A specification cannot have zero performance inputs.'), {obj.specs.in}); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'specs.in can only contain numeric values.'), {obj.specs.in}); 
            cellfun(@(x) obj.watchdog.assertCond(all(mod(x,1)==0) && all(x>0), 'specs.in can only contain integer values.'), {obj.specs.in});
            cellfun(@(x) obj.watchdog.assertNonEmpty(x, 'A specification cannot have zero performance outputs.'), {obj.specs.out}); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'specs.out can only contain numeric values.'), {obj.specs.out});
            cellfun(@(x) obj.watchdog.assertCond(all(mod(x,1)==0) && all(x>0), 'specs.out can only contain integer values.'), {obj.specs.out});
            cellfun(@(x) obj.watchdog.assertNonEmpty(x,'A specification requires a weight. Set to 0 to handle it as a normalized constraint.'), {obj.specs.weight}); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'specs.weight can only contain numeric values.'), {obj.specs.weight}); 
            cellfun(@(x) obj.watchdog.assertCond(all(x>=0) && isscalar(x), 'specs.weight can only positive scalars or 0.'), {obj.specs.weight});
            obj.watchdog.assertField(obj.region, 'M', 'region requires a field ''M''.'); 
            obj.watchdog.assertField(obj.region, 'L', 'region requires a field ''L''.'); 
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'region.M can only contain numeric values.'), {obj.region.M});
            cellfun(@(x) obj.watchdog.assertType(x, 'double', 'region.L can only contain numeric values.'), {obj.region.L});
            cellfun(@(x) obj.watchdog.assertCond(isreal(x) && issymmetric(x), 'region.L can only contain real symmetric matrices.'), {obj.region.L});
            cellfun(@(x) obj.watchdog.assertCond(isreal(x), 'region.M can only contain real matrices.'), {obj.region.M}); 
            cellfun(@(x,y) obj.watchdog.assertCond(size(x,1)==size(y,1) && issymmetric(y), 'region.M and region.L must have the same number of rows. In addition, L must be symmetric.'), {obj.region.M}, {obj.region.L});
            
            obj.watchdog.assertCond(max(vertcat(obj.specs.in))<=size(obj.Wi,2), 'A performance input exceeding the inputs of the input weighting filter was requested.'); 
            obj.watchdog.assertCond(max(vertcat(obj.specs.out))<=size(obj.Wo,1), 'A performance output exceeding the outputs of the output weighting filter was requested.'); 
            obj.watchdog.warnCond(all(ismember(1:size(obj.Wi,2),vertcat(obj.specs.in))),'Some performance inputs do not appear in the specifications. They will therefore be ignored.');
            obj.watchdog.warnCond(all(ismember(1:size(obj.Wo,1),vertcat(obj.specs.out))),'Some performance outputs do not appear in the specifications. They will therefore be ignored.');
        end
        
        function obj = sortchannels(obj)
        % SORTCHANNELS Sorts the specifications and its performance inputs and outputs
            obj.Wi = obj.Wi(:,vertcat(obj.specs.in));
            obj.Wo = obj.Wo(vertcat(obj.specs.out),:); 
            chi = 0; cho = 0;
            for k=1:length(obj.specs)
                obj.specs(k).in = (chi+1):(chi+length(obj.specs(k).in)); chi = chi+length(obj.specs(k).in);
                obj.specs(k).out = (cho+1):(cho+length(obj.specs(k).out)); cho = cho+length(obj.specs(k).out); 
            end
        end
        
        function checkassumptions(obj)
        % CHECKASSUMPTIONS Check the assumptions for the extended H-infinity control problem

            % conditions on P
            obj.watchdog.assertCond(sshinfcd.util.isdtcb(obj.P) && sshinfcd.util.isstbz(obj.P), 'P should be finite dynamics detectable and stabilizable.'); 
            
            % conditions on Wi 
            obj.watchdog.assertCond(sshinfcd.util.isctrb(obj.Wi), 'Wi should be controllable.');
            
            % conditions on Wo
            obj.watchdog.assertCond(sshinfcd.util.isobsv(obj.Wo), 'Wo should be observable.');
            
            % interdependent conditions
            WoP = blkdiag(obj.Wo,eye(obj.ny()))*obj.P;
            PWi = obj.P*blkdiag(obj.Wi,eye(obj.nu()));
            WoP = obj.watchdog.executeSilently('minreal(WoP)','An error occured while calculating a minimal realization of the combination of P and Wo (Phat).'); 
            PWi = obj.watchdog.executeSilently('minreal(PWi)','An error occured while calculating a minimal realization of P and Wi (Ptilde).'); 
            obj.watchdog.assertCond(sshinfcd.util.isstbz(WoP(:,(obj.nv())+1:end)),'The combination of P and Wo (Phat) should be finite dynamics stabilizable.');
            obj.watchdog.assertCond(sshinfcd.util.isdtcb(PWi((obj.ne()+1):end,:)),'The combination of P and Wi (Ptilde) should be finite dynamics detectable.');
            obj.watchdog.assertCond(isempty(intersect(pole(obj.Wi),pole(obj.Wo))),'Wi and Wo cannot share poles.'); 
        end
        
        function obj = standardfilters(obj)
        % STANDARDFILTERS Incorporates stable finite dynamic modes from Wi and Wo into P
            % input weight
            [WIS,WIUS] = obj.watchdog.executeSilently('sshinfcd.util.stabsep(obj.Wi,obj.region)','An error occured while separating the D-stable and the D-unstable modes from the input filter. Check the logbook for more details.');
            
            % output weight
            [WOS,WOUS] = obj.watchdog.executeSilently('sshinfcd.util.stabsep(obj.Wo,obj.region)','An error occured while separating the D-stable and the D-unstable modes from the output filter. Check the logbook for more details.');

            % merge with plant
            obj.P = blkdiag([WOS ; eye(obj.ne())],eye(obj.ny()))*obj.P*blkdiag([eye(obj.nv()) WIS],eye(obj.nu())); 
            obj.Wi = [WIUS ; eye(obj.nw())]; 
            obj.Wo = [eye(obj.nz()) WOUS]; 
        end
        
        %% Different canonical realizations
        
        function obj = minreal(obj)
        % MINREAL Returns a minimal realization of P, Wi and Wo
            obj.P = obj.watchdog.executeSilently('minreal(obj.P)','An error occured while calculating the minimal realization of the plant P. Check the logbook for more details.');
            obj.Wi = obj.watchdog.executeSilently('minreal(obj.Wi)','An error occured while calculating the minimal realization of the input filter Wi. Check the logbook for more details.');
            obj.Wo = obj.watchdog.executeSilently('minreal(obj.Wo)','An error occured while calculating the minimal realization of the output filter Wo. Check the logbook for more details.');
        end
        
        function obj = balreal(obj)
        % BALREAL Balances the stable finite dynamics of P
            obj.P = balreal(obj.P);
        end
        
        %% Sampling time
        
        function Ts = Ts(obj)
        % TS Returns the sampling time of the plant and the filters
            Ts = obj.P.Ts;
        end
        
        %% Dimensions
        
        function nv = nv(obj)
        % NV Returns number of performance inputs of P = number of outputs of Wi
            nv = size(obj.Wi,1);
        end
        
        function nw = nw(obj)
        % NW Returns number of inputs of Wi
            nw = size(obj.Wi,2);
        end
        
        function ne = ne(obj)
        % NE Returns number of inputs of Wo = number of performance outputs of P
            ne = size(obj.Wo,2);
        end
        
        function nz = nz(obj)
        % NZ Returns number of outputs of Wo
            nz = size(obj.Wo,1);
        end
        
        function nu = nu(obj)
        % NU Returns number of control inputs
            nu = size(obj.P,2)-obj.nv(); 
        end
        
        function ny = ny(obj)
        % NY Returns number of measured outputs
            ny = size(obj.P,1)-obj.ne(); 
        end
        
        function ni = ni(obj)
        % NI Returns order of the input filter
            ni = size(obj.Ai,1);
        end
        
        function no = no(obj)
        % NO Returns order of the output filter
            no = size(obj.Ao,1);
        end
        
        function n = n(obj)
        % N Returns order of the plant
            n = size(obj.A,1);
        end
        
        %% Matrices of P        
        
        function A = A(obj)
        % A Returns the A-matrix of the generalized plant
            A = obj.P.A;
        end
        
        function Bv = Bv(obj)
        % BV Returns the Bv-matrix of the generalized plant
            Bv = obj.P.B(:,1:obj.nv());
        end
        
        function Bu = Bu(obj)
        % BU Returns the Bv-matrix of the generalized plant
            Bu = obj.P.B(:,(obj.nv()+1):end);
        end
        
        function Ce = Ce(obj)
        % CE Returns the Ce-matrix of the generalized plant
            Ce = obj.P.C(1:obj.ne(),:);
        end
        
        function Cy = Cy(obj)
        % CY Returns the Cy-matrix of the generalized plant
            Cy = obj.P.C((obj.ne()+1):end,:);
        end
        
        function Dev = Dev(obj)
        % DEV Returns the Dev-matrix of the generalized plant
            Dev = obj.P.D(1:obj.ne(),1:obj.nv());
        end
        
        function Deu = Deu(obj)
        % DEU Returns the Deu-matrix of the generalized plant
            Deu = obj.P.D(1:obj.ne(),(obj.nv()+1):end);
        end
        
        function Dyv = Dyv(obj)
        % DYV Returns the Dyv-matrix of the generalized plant
            Dyv = obj.P.D((obj.ne()+1):end,1:obj.nv());
        end
        
        function Dyu = Dyu(obj)
        % DYU Returns the Dyu-matrix of the generalized plant
            Dyu = obj.P.D((obj.ne()+1):end,(obj.nv()+1):end);
        end
        
        %% Matrices of Wi
        
        function Ai = Ai(obj)
        % AI Returns the A-matrix of the input filter
            Ai = obj.Wi.A;
        end
            
        function Bi = Bi(obj)
        % BI Returns the B-matrix of the input filter
            Bi = obj.Wi.B;
        end
        
        function Ci = Ci(obj)
        % CI Returns the C-matrix of the input filter
            Ci = obj.Wi.C;
        end
        
        function Di = Di(obj)
        % DI Returns the D-matrix of the input filter
            Di = obj.Wi.D;
        end
        
        %% Matrices of Wo
        
        function Ao = Ao(obj)
        % AO Returns the A-matrix of the output filter
            Ao = obj.Wo.A;
        end
        
        function Bo = Bo(obj)
        % BO Returns the B-matrix of the output filter
            Bo = obj.Wo.B;
        end
        
        function Co = Co(obj)
        % CO Returns the C-matrix of the output filter
            Co = obj.Wo.C;
        end
        
        function Do = Do(obj)
        % DO Returns the D-matrix of the output filter
            Do = obj.Wo.D;
        end
        
        %% Matrices of Phat
        
        function Ahat = Ahat(obj)
        % EHAT Returns the E-matrix of Phat
            Ahat = [obj.Ao()                                  obj.Bo()*obj.Ce(); 
                    zeros(size(obj.A(),1),size(obj.Ao(),2))   obj.A()          ];
        end
        
        function Bvhat = Bvhat(obj)
        % BVHAT Returns the Bv-matrix of Phat
            Bvhat = [obj.Bo()*obj.Dev();
                     obj.Bv()          ];
        end
        
        function Buhat = Buhat(obj)
        % BUHAT Returns the Bu-matrix of Phat
            Buhat = [obj.Bo()*obj.Deu();
                     obj.Bu()          ];
        end
        
        function Czhat = Czhat(obj)
        % CZHAT Returns the Cz-matrix of Phat
            Czhat = [obj.Co()  obj.Do()*obj.Ce()];
        end
        
        function Cyhat = Cyhat(obj)
        % CYHAT Returns the Cy-matrix of Phat
            Cyhat = [zeros(size(obj.Cy(),1),size(obj.Ao(),2)) obj.Cy()];
        end
        
        function Dzvhat = Dzvhat(obj)
        % DZVHAT Returns the Dzv-matrix of Phat
            Dzvhat = obj.Do()*obj.Dev();
        end
        
        function Dzuhat = Dzuhat(obj)
        % DZUHAT Returns the Dzu-matrix of Phat
            Dzuhat = obj.Do()*obj.Deu();
        end
        
        function Dyvhat = Dyvhat(obj)
        % DYVHAT Returns the Dyv-matrix of Phat
            Dyvhat = obj.Dyv(); 
        end
        
        function Phat = Phat(obj)
        % PHAT Returns the DSS realization of Phat 
            Phat = ss(obj.Ahat(), [obj.Bvhat() obj.Buhat()], [obj.Czhat() ; obj.Cyhat()], [obj.Dzvhat() obj.Dzuhat() ; obj.Dyvhat() obj.Dyu()], obj.Ts());
        end
        
        %% Matrices of Ptilde
        
        function Atilde = Atilde(obj)
        % ATILDE Returns the A-matrix of Ptilde
            Atilde = [obj.A()                                   obj.Bv()*obj.Ci();
                      zeros(size(obj.Ai(),1),size(obj.A(),2))   obj.Ai()         ];
        end
        
        function Bwtilde = Bwtilde(obj)
        % BWTILDE Returns the Bw-matrix of Ptilde
            Bwtilde = [obj.Bv()*obj.Di();
                       obj.Bi()         ];
        end
        
        function Butilde = Butilde(obj)
        % BUTILDE Returns the Bu-matrix of Ptilde
            Butilde = [obj.Bu()                                ; 
                       zeros(size(obj.Ai(),2),size(obj.Bu(),2))];
        end
        
        function Cetilde = Cetilde(obj)
        % CETILDE Returns the Ce-matrix of Ptilde
            Cetilde = [obj.Ce()  obj.Dev()*obj.Ci()];
        end
        
        function Cytilde = Cytilde(obj)
        % CYTILDE Returns the Cy-matrix of Ptilde
            Cytilde = [obj.Cy()  obj.Dyv()*obj.Ci()];
        end
        
        function Dewtilde = Dewtilde(obj)
        % DEWTILDE Returns the Dew-matrix of Ptilde
            Dewtilde = obj.Dev()*obj.Di(); 
        end
        
        function Deutilde = Deutilde(obj)
        % DEUTILDE Returns the Deu-matrix of Ptilde
            Deutilde = obj.Deu(); 
        end
        
        function Dywtilde = Dywtilde(obj)
        % DYWTILDE Returns the Deu-matrix of Ptilde
            Dywtilde = obj.Dyv()*obj.Di(); 
        end
        
        function Ptilde = Ptilde(obj)
        % PTILDE Returns the DSS realization of Ptilde 
            Ptilde = ss(obj.Atilde(), [obj.Bwtilde() obj.Butilde()], [obj.Cetilde() ; obj.Cytilde()], [obj.Dewtilde() obj.Deutilde() ; obj.Dywtilde() obj.Dyu()], obj.Ts());
        end
        
        %% Matrices of Pbar
        
        function Abar = Abar(obj)
        % ABAR Returns the A-matrix of Pbar
            Abar = [obj.Ao()                                      obj.Bo()*obj.Cetilde(); 
                    zeros(size(obj.Atilde(),1),size(obj.Ao(),2))  obj.Atilde()          ];
        end
        
        function Bwbar = Bwbar(obj)
        % BWBAR Returns the Bw-matrix of Pbar
            Bwbar = [obj.Bo()*obj.Dewtilde(); 
                     obj.Bwtilde()          ];
        end
        
        function Bubar = Bubar(obj)
        % BUBAR Returns the Bu-matrix of Pbar
            Bubar = [obj.Bo()*obj.Deutilde(); 
                     obj.Butilde()          ]; 
        end
        
        function Czbar = Czbar(obj)
        % CZBAR Returns the Cz-matrix of Pbar
            Czbar = [obj.Co() obj.Do()*obj.Cetilde()]; 
        end
        
        function Cybar = Cybar(obj)
        % CYBAR Returns the Cy-matrix of Pbar
            Cybar = [zeros(size(obj.Cytilde(),1),size(obj.Ao(),2)) obj.Cytilde()]; 
        end
        
        function Dzwbar = Dzwbar(obj)
        % DZWBAR Returns the Dzw-matrix of Pbar
            Dzwbar = obj.Do()*obj.Dewtilde(); 
        end
        
        function Dzubar = Dzubar(obj)
        % DZUBAR Returns the Dzu-matrix of Pbar
            Dzubar = obj.Do()*obj.Deutilde(); 
        end
        
        function Dywbar = Dywbar(obj)
        % DYWBAR Returns the Dyw-matrix of Pbar
            Dywbar = obj.Dywtilde(); 
        end
        
        function Pbar = Pbar(obj)
        % PBAR Returns the DSS realization of Pbar 
            Pbar = ss(obj.Abar(), [obj.Bwbar() obj.Bubar()], [obj.Czbar() ; obj.Cybar()], [obj.Dzwbar() obj.Dzubar() ; obj.Dywbar() obj.Dyu()], obj.Ts());
        end
        
        %% Elimination of the D-unstable modes in the weighting filters
        
        function [GAMMAi,PIihat,Zi] = ifseparator(obj)
        % IFSEPARATOR Returns the separator variables for the input filter
        % Returns the algebraic variables GAMMAi, PIihat and Zi that
        % eliminate the input weighting filter dynamics from the
        % generalized plant and are required for the reconstruction of the
        % controller that renders the closed-loop input-output stable.
        
            % solve Sylvester equation
            A = [obj.Ahat() obj.Buhat() ; obj.Czhat() obj.Dzuhat()];
            B = -obj.Ai();
            C = [obj.Bvhat*obj.Ci() ; obj.Dzvhat()*obj.Ci()]; 
            X = obj.watchdog.executeSilently('sshinfcd.util.strucsylvester(A,B,C,obj.n()+obj.no());', 'Elimination of the input filter dynamics failed. Check the logbook for more details. Do your performance outputs depend on the control inputs?');
            
            % extract variables
            PIihat = X(1:(obj.no()+obj.n()),:);
            GAMMAi = X((obj.n()+obj.no()+1):end,:);
            Zi = obj.Cy()*PIihat((obj.no()+1):end,:) - obj.Dyv()*obj.Ci();
        end
        
        function [GAMMAo,PIotilde,Zotilde] = ofseparator(obj)
        % IFSEPARATOR Returns the separator variables for the input filter
        % Returns the algebraic variables GAMMAo, PIotilde and Zotilde 
        % that eliminate the output weighting filter dynamics from the
        % generalized plant and are required for the reconstruction of the
        % controller that renders the closed-loop input-output admissible.
        
            % solve Sylvester equation
            A = [obj.Atilde() obj.Bwtilde() ; obj.Cytilde() obj.Dywtilde()];
            B = -obj.Ao();
            C = [obj.Bo()*obj.Cetilde() obj.Bo()*obj.Dewtilde()];
            X = obj.watchdog.executeSilently('sshinfcd.util.strucsylvester(transpose(A),transpose(B),transpose(C),obj.ni()+obj.n());', 'Elimination of the output filter dynamics failed. Check the logbook for more details. Do your performance inputs affect the measurement outputs?');
            X = X'; 
            
            % get separation variables of the input filter
            [~,PIihat,~] = obj.ifseparator(); 
            
            % extract variables
            PIotilde = X(:,1:obj.n()+obj.ni());
            GAMMAo = X(:,(obj.n()+obj.ni()+1):end);
            Zotilde = [PIotilde(:,1:obj.n())*obj.Bu() - obj.Bo()*obj.Deu(), PIihat(1:obj.no(),:) + PIotilde(:,(obj.n()+1):end) - PIotilde(:,1:obj.n())*PIihat((obj.no()+1):end,:)];
        end
        
        function Buiobar = Buiobar(obj)
        % BUIOBAR Returns the Bu-matrix of the extended plant after elimination of the filter dynamics (Piobar)
            [~,PIihat,~] = ifseparator(obj);
            Buiobar = [-PIihat        obj.Buhat()                         ;
                       eye(obj.ni())  zeros(obj.ni(), size(obj.Buhat(),2))];
        end
        
        function Dzuiobar = Dzuiobar(obj)
        % DZUIOBAR Returns the Dzu-matrix of the extended plant after elimination of the filter dynamics (Piobar)
            Dzuiobar = [zeros(size(obj.Dzubar(),1), obj.ni()) obj.Dzubar()];
        end
        
        function Cyiobar = Cyiobar(obj)
        % CYIOBAR Returns the Cy-matrix of the extended plant after elimination of the filter dynamics (Piobar)
            [~,PIotilde,~] = ofseparator(obj);
            Cyiobar = [eye(obj.no())                           -PIotilde;
                       zeros(size(obj.Cytilde(),1), obj.no())  obj.Cytilde()];
        end
        
        function Dywiobar = Dywiobar(obj)
        % DYWIOBAR Returns the Dyw-matrix of the extended plant after elimination of the filter dynamics (Piobar)
            Dywiobar = [zeros(obj.no(), size(obj.Dywbar(),2)) ; obj.Dywbar()];
        end
        
        function Piobar = Piobar(obj)
        % PIOBAR Returns the DSS realization of the extended plant after elimination of the filter dynamics (Piobar)
            obj.watchdog.warnCond(0,'Feedthrough part Dyu() is set to zero to construct Piobar.');
            Piobar = ss(obj.Abar(), [obj.Bwbar() obj.Buiobar()], [obj.Czbar() ; obj.Cyiobar()], [obj.Dzwbar() obj.Dzuiobar() ; obj.Dywiobar() zeros(obj.ny()+obj.no(),obj.nu()+obj.ni())], obj.Ts());
        end
       
        %% Helper variables for the LMI formulation
        
        function XI = XI(obj)
        % XI Returns the off-diagonal stability LMI blocks
            [~,PIotilde,~] = obj.ofseparator();
            [~,PIihat,~] = obj.ifseparator(); 
            XI = [PIotilde(:,1:obj.n())             PIihat(1:obj.no(),:) + PIotilde(:,(obj.n()+1):end);
                  eye(obj.n())                      PIihat((obj.no()+1):end,:)                        ];
        end
        
        function PIi = PIi(obj)
        % PIi Returns PIi 
            [~,PIihat,~] = obj.ifseparator(); 
            PIi = PIihat((obj.no()+1):end,:);
        end
        
        function PIo = PIo(obj)
        % PIo Returns PIo 
            [~,PIotilde,~] = obj.ofseparator();
            PIo = PIotilde(:,1:obj.n());
        end
        
        function Czlmi = Czlmi(obj)
        % CZLMI Returns the Cz block of the performance LMI
            [~,PIotilde,~] = obj.ofseparator();
            Czlmi = obj.Co()*PIotilde + obj.Do()*obj.Cetilde(); 
        end
        
        function Bwlmi = Bwlmi(obj)
        % BWLMI Returns the Bw block of the performance LMI 
            [~,PIihat,~] = obj.ifseparator();
            Bwlmi = PIihat*obj.Bi() + obj.Bvhat()*obj.Di(); 
        end
        
        function PoTT = PoTTlmi(obj)
        % POTTLMI Returns the Po'*T block of the performance LMI
            [~,PIihat,~] = obj.ifseparator();
            PoTT = [eye(obj.n()+obj.no()) PIihat];
        end
        
        function SPi = SPilmi(obj)
        % SPILMI Returns the S*Pi block of the performance LMI
            [~,PIotilde,~] = obj.ofseparator();
            SPi = [PIotilde ; eye(obj.n()+obj.ni())];
        end
        
        function PoTTASPi = PoTTASPilmi(obj)
        % PoTTASPi Returns the Po'*T*A*S*Pi block of the performance LMI
            PoTTASPi = obj.PoTTlmi()*obj.Abar()*obj.SPilmi(); 
        end
        
    end
    
end