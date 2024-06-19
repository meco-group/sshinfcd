function [K, gamma] = nlcov(prob,opts)
% NLCOV Solves the LMI formulation of the standard problem in state-space form
%
% [K,gamma] = NLCOV(prob) returns a descriptor realization K of the
% controller that solves the SDP formulation of the standard problem,
% i.e. for the plant Pio (after elimination of impulsive and unstable 
% weighting filters)

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

    assert(prob.Ts==0,'Discrete-time problems not supported yet.');
    
    % 1. Solve the synthesis problem
    switch lower(opts.synthesis.gammasolver)
        case 'yalmip'
            sdpvars = ct_synthesis_YALMIP(prob,opts);
        case 'cvx'
            sdpvars = ct_synthesis_CVX(prob,opts);
        case 'lmilab'
            sdpvars = ct_synthesis_LMILAB(prob,opts);
        otherwise
            prob.watchdog.error('LMI parser for synthesis problem not defined.'); 
    end
    
    % 2. Reconstruct the controller 
    K = ct_controller_reconstruction(sdpvars,prob);
    
    % 3. Return gamma or gamma^2 (depending on Lyapunov shaping formulation)
    gamma = sdpvars.sqgamma; 

end

%% CONTINUOUS TIME - SYNTHESIS

function sdpvars = ct_synthesis_YALMIP(prob,opts)
    
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % define main variables
        % Lyapunov: X, Y
        % controller: Kc, Mc, Nc, Dc
        X = sdpvar(prob.n()+prob.ni(),prob.n()+prob.ni(),'symmetric');
        Y = sdpvar(prob.n()+prob.no(),prob.n()+prob.no(),'symmetric');
        Kc = sdpvar(prob.n()+prob.ni(),prob.n()+prob.no(),'full');
        Mc = sdpvar(prob.n()+prob.ni(),prob.ny(),'full');
        Nc = sdpvar(prob.nu(),prob.n()+prob.no(),'full');
        Dc = sdpvar(prob.nu(),prob.ny(),'full'); 
        
        % gamma
        sqgamma = sdpvar(nobj,1,'full');
        
    % construct matrices with slack variables
        Gammaw = sdpvar(prob.nw(),prob.nw(),'symmetric');
        Gammaz = sdpvar(prob.nz(),prob.nz(),'symmetric');
        Dzw = sdpvar(prob.nz(),prob.nw(),'full');
        Dzwcst = prob.Dzwbar(); 
        Dzucst = prob.Dzubar();
        Dywcst = prob.Dywbar(); 
        
        % set constant parts & replace optimization variables
        j = 1;
        for i=1:length(prob.specs)
        	if prob.specs(i).weight==0
                Gammaz(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out)); 
                Gammaw(prob.specs(i).in,prob.specs(i).in) = eye(length(prob.specs(i).in)); 
            else
                switch opts.synthesis.lyapunovshape 
                    case 1
                        Gammaz(prob.specs(i).out,prob.specs(i).out) = sqgamma(j)*eye(length(prob.specs(i).out));
                        Gammaw(prob.specs(i).in,prob.specs(i).in) = sqgamma(j)*eye(length(prob.specs(i).in));
                        j = j+1;
                    case 2
                        Gammaz(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out));
                        Gammaw(prob.specs(i).in,prob.specs(i).in) = sqgamma(j)*eye(length(prob.specs(i).in));
                        j = j+1;
                    case 3
                        Gammaz(prob.specs(i).out,prob.specs(i).out) = sqgamma(j)*eye(length(prob.specs(i).out));
                        Gammaw(prob.specs(i).in,prob.specs(i).in) = eye(length(prob.specs(i).in));
                        j = j+1;
                    otherwise
                        prob.watchdog.error('The Lyapunov shaping format can only be of type 1 (-gamma,-gamma), type 2 (-gamma^2,-1) or type 3 (-1,-gamma^2).'); 
                end
            end
            Dzw(prob.specs(i).out,prob.specs(i).in) = Dzwcst(prob.specs(i).out,prob.specs(i).in) + Dzucst(prob.specs(i).out,:)*Dc*Dywcst(:,prob.specs(i).in);
        end
        
    % formulate LMIs
        % STABILITY 
        STAB = [X          prob.XI()'; 
                prob.XI()  Y         ];

        % H-INFINITY PERFORMANCE 
        PERF =     [prob.Ahat()*Y+Y*prob.Ahat()'+prob.Buhat()*Nc+Nc'*prob.Buhat()',    prob.PoTTASPilmi()+prob.Buhat()*Dc*prob.Cytilde()+Kc',                      prob.Bwlmi()+prob.Buhat()*Dc*prob.Dywtilde(),        Y*prob.Czhat()'+Nc'*prob.Dzuhat()'              ;
                    prob.PoTTASPilmi()'+Kc+prob.Cytilde()'*Dc'*prob.Buhat()',          X*prob.Atilde()+prob.Atilde()'*X+Mc*prob.Cytilde()+prob.Cytilde()'*Mc',     X*prob.Bwtilde()+Mc*prob.Dywtilde(),                 prob.Czlmi()'+prob.Cytilde()'*Dc'*prob.Dzuhat()';
                    prob.Bwlmi()'+prob.Dywtilde()'*Dc'*prob.Buhat()',                  prob.Bwtilde()'*X+prob.Dywtilde()'*Mc',                                     -Gammaw,                                             Dzw'                                            ;
                    prob.Czhat()*Y+prob.Dzuhat()*Nc,                                   prob.Czlmi()+prob.Dzuhat()*Dc*prob.Cytilde(),                               Dzw,                                                 -Gammaz                                         ]; 
 
        % D-STABILITY
        DSTAB = {};
        for i=1:length(prob.region)
            if ~(isempty(prob.region(i).M) && isempty(prob.region(i).L))
                DSTAB{i} = kron(prob.region(i).L, STAB) + 2*kron(prob.region(i).M,[prob.Ahat()*Y+prob.Buhat()*Nc,     prob.PoTTASPilmi()+prob.Buhat()*Dc*prob.Cytilde();
                                                                                   Kc,                                X*prob.Atilde()+Mc*prob.Cytilde()                ]);
            end
        end
                     
   % optimize
   obj = vertcat(prob.specs(isobj).weight)'*sqgamma;
   cstr = [0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0;
           0.5*(PERF+PERF') - opts.synthesis.zerotol*eye(size(PERF)) <= 0];
           for i=1:length(DSTAB)
               cstr = [cstr; 
                       0.5*(DSTAB{i}+DSTAB{i}') - opts.synthesis.zerotol*eye(size(DSTAB{i})) <= 0];
           end
   optimize(cstr,obj,opts.yalmip); 
   
   % put into the output structure
   sdpvars.X = double(X); 
   sdpvars.Y = double(Y); 
   sdpvars.Kc = double(Kc); 
   sdpvars.Mc = double(Mc);
   sdpvars.Nc = double(Nc);
   sdpvars.Dc = double(Dc); 
   sdpvars.sqgamma = double(sqgamma); 
   
end

function sdpvars = ct_synthesis_CVX(prob, opts)
   
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
        
    % start CVX environment    
    cvx_begin SDP
        
        % Lyapunov: X, Y
        % controller: Kc, Mc, Nc, Dc
        variable X(prob.n()+prob.ni(),prob.n()+prob.ni()) symmetric
        variable Y(prob.n()+prob.no(),prob.n()+prob.no()) symmetric
        variable Kc(prob.n()+prob.ni(),prob.n()+prob.no()) 
        variable Mc(prob.n()+prob.ni(),prob.ny()) 
        variable Nc(prob.nu(),prob.n()+prob.no()) 
        variable Dc(prob.nu(),prob.ny()) 
        
        % gamma
        variable sqgamma(nobj,1)

        % construct structured matrix variables 
        j = 1;
        if prob.specs(1).weight==0
            Gammaw = eye(length(prob.specs(1).in));
            Gammaz = eye(length(prob.specs(1).out));
        else
            switch opts.synthesis.lyapunovshape 
                case 1
                    Gammaw = sqgamma(j)*eye(length(prob.specs(1).in));
                    Gammaz = sqgamma(j)*eye(length(prob.specs(1).out));
                case 2
                    Gammaw = sqgamma(j)*eye(length(prob.specs(1).in));
                    Gammaz = eye(length(prob.specs(1).out));
                case 3
                    Gammaw = eye(length(prob.specs(1).in));
                    Gammaz = sqgamma(j)*eye(length(prob.specs(1).out));
                otherwise
                    prob.watchdog.error('The Lyapunov shaping format can only be of type 1 (-gamma,-gamma), type 2 (-gamma^2,-1) or type 3 (-1,-gamma^2).'); 
            end
            j = j+1;
        end
        Dzwcst = prob.Dzwbar(); 
        Dzucst = prob.Dzubar();
        Dywcst = prob.Dywbar(); 
        Dzw = Dzwcst(prob.specs(1).out,prob.specs(1).in);
        for i=2:length(prob.specs)
            eval(['variable slackw' num2str(i) '(size(Gammaw,1),length(prob.specs(i).in))']);
            eval(['variable slackz' num2str(i) '(size(Gammaz,1),length(prob.specs(i).out))']);
            
            if prob.specs(i).weight==0 
                eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' eye(length(prob.specs(i).out))];']);
                eval(['Gammaw = [Gammaw slackw' num2str(i) '; slackw' num2str(i) ''' eye(length(prob.specs(i).in))];']);
            else
                switch opts.synthesis.lyapunovshape
                    case 1
                        eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' sqgamma(j)*eye(length(prob.specs(i).out))];']);
                        eval(['Gammaw = [Gammaw slackw' num2str(i) '; slackw' num2str(i) ''' sqgamma(j)*eye(length(prob.specs(i).in))];']);
                    case 2
                        eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' sqgamma(j)*eye(length(prob.specs(i).out))];']);
                        eval(['Gammaw = [Gammaw slackw' num2str(i) '; slackw' num2str(i) ''' eye(length(prob.specs(i).in))];']);
                    case 3
                        eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' eye(length(prob.specs(i).out))];']);
                        eval(['Gammaw = [Gammaw slackw' num2str(i) '; slackw' num2str(i) ''' sqgamma(j)*eye(length(prob.specs(i).in))];']);
                    otherwise
                        prob.watchdog.error('The Lyapunov shaping format can only be of type 1 (-gamma,-gamma), type 2 (-gamma^2,-1) or type 3 (-1,-gamma^2).'); 
                end
                j = j+1;
             end
            
            % Dzw
            eval(['variable slackdh' num2str(i) '(size(Dzw,1),length(prob.specs(i).in))']);
            eval(['variable slackdv' num2str(i) '(length(prob.specs(i).out),size(Dzw,2))']);
            eval(['Dzw = [Dzw slackdh' num2str(i) '; slackdv' num2str(i) ' (Dzwcst(prob.specs(i).out,prob.specs(i).in)+Dzucst(prob.specs(i).out,:)*Dc*Dywcst(:,prob.specs(i).in))];']);
        end        
        
        % formulate LMIs
            % STABILITY 
            STAB = [X          prob.XI()'; 
                    prob.XI()  Y         ];

            % H-INFINITY PERFORMANCE 
            PERF =     [prob.Ahat()*Y+Y*prob.Ahat()'+prob.Buhat()*Nc+Nc'*prob.Buhat()',    prob.PoTTASPilmi()+prob.Buhat()*Dc*prob.Cytilde()+Kc',                      prob.Bwlmi()+prob.Buhat()*Dc*prob.Dywtilde(),        Y*prob.Czhat()'+Nc'*prob.Dzuhat()'              ;
                        prob.PoTTASPilmi()'+Kc+prob.Cytilde()'*Dc'*prob.Buhat()',          X*prob.Atilde()+prob.Atilde()'*X+Mc*prob.Cytilde()+prob.Cytilde()'*Mc',     X*prob.Bwtilde()+Mc*prob.Dywtilde(),                 prob.Czlmi()'+prob.Cytilde()'*Dc'*prob.Dzuhat()';
                        prob.Bwlmi()'+prob.Dywtilde()'*Dc'*prob.Buhat()',                  prob.Bwtilde()'*X+prob.Dywtilde()'*Mc',                                     -Gammaw,                                             Dzw'                                            ;
                        prob.Czhat()*Y+prob.Dzuhat()*Nc,                                   prob.Czlmi()+prob.Dzuhat()*Dc*prob.Cytilde(),                               Dzw,                                                 -Gammaz                                         ]; 

            % D-STABILITY
            DSTAB = {};
            for i=1:length(prob.region)
                if ~(isempty(prob.region(i).M) && isempty(prob.region(i).L))
                    DSTAB{i} = kron(prob.region(i).L, STAB) + 2*kron(prob.region(i).M,[prob.Ahat()*Y+prob.Buhat()*Nc,     prob.PoTTASPilmi()+prob.Buhat()*Dc*prob.Cytilde();
                                                                                       Kc,                                X*prob.Atilde()+Mc*prob.Cytilde()                ]);
                    DSTAB{i} = (0.5*(DSTAB{i}+DSTAB{i}') - opts.synthesis.zerotol*eye(size(DSTAB{i})) <= 0);
                end
            end
                
        % optimize
        minimize([prob.specs(isobj).weight]*sqgamma)
        subject to
            0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0
            0.5*(PERF+PERF') - opts.synthesis.zerotol*eye(size(PERF)) <= 0
            DSTAB{:}
            
         eval(['cvx_solver ' opts.cvx.solver]);
         eval(['cvx_precision ' opts.cvx.precision]);
         eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end
    
   % put into the output structure
   sdpvars.X = X;
   sdpvars.Y = Y; 
   sdpvars.Kc = Kc; 
   sdpvars.Mc = Mc;
   sdpvars.Nc = Nc;
   sdpvars.Dc = Dc; 
   sdpvars.sqgamma = sqgamma; 
    
end

function sdpvars = ct_synthesis_LMILAB(prob,opts)
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % define main variables
    setlmis([]);
    
        % Lyapunov: X, Y
        % controller: Kc, Mc, Nc, Dc
        X = lmivar(1,[prob.n()+prob.ni(),1]);
        Y = lmivar(1,[prob.n()+prob.no(),1]); 
        Kc = lmivar(2,[prob.n()+prob.ni(),prob.n()+prob.no()]); 
        Mc = lmivar(2,[prob.n()+prob.ni(),prob.ny()]); 
        Nc = lmivar(2,[prob.nu(),prob.n()+prob.no()]); 
        Dc = lmivar(2,[prob.nu(),prob.ny()]); 
        
        % gamma
        [sqgamma,nvar,ssqgamma] = lmivar(2,[nobj,1]);
    
    % construct matrices with slack variables
        % helper function
        function struc = offdiagonal(Nz,Nw,full)
            struc = zeros(sum(Nz),sum(Nw));
            Nw = [0 cumsum(Nw)]; Nz = [0 cumsum(Nz)]; 
            for ii=1:length(Nw)-1
                struc((Nz(ii+1)+1):end,(Nw(ii)+1):Nw(ii+1)) = 1;
                if nargin>=3 && full
                    struc(1:Nz(ii),(Nw(ii)+1):Nw(ii+1)) = 1;
                end
            end
        end

        % matrix structures
        iw = find(offdiagonal(cellfun(@length, {prob.specs.in}),cellfun(@length, {prob.specs.in}))); 
        iz = find(offdiagonal(cellfun(@length, {prob.specs.out}),cellfun(@length, {prob.specs.out}))); 
        id = find(offdiagonal(cellfun(@length, {prob.specs.out}),cellfun(@length, {prob.specs.in}),1)); 
        
        % create slack variables 
        sGammaw = zeros(prob.nw()); sGammaw(iw) = (1:length(iw))+nvar; sGammaw = sGammaw+sGammaw';
        sGammaz = zeros(prob.nz()); sGammaz(iz) = (1:length(iz))+nvar+length(iw); sGammaz = sGammaz+sGammaz';
        sD = zeros(prob.nz(),prob.nw()); sD(id) = (1:length(id))+nvar+length(iw)+length(iz); 
        
        % append optimization variables
        j = 1;
        Gammawcst = zeros(prob.nw()); 
        Gammazcst = zeros(prob.nz()); 
        for i=1:length(prob.specs)
        	if prob.specs(i).weight==0
                Gammawcst(prob.specs(i).in,prob.specs(i).in) = eye(length(prob.specs(i).in)); 
                Gammazcst(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out)); 
            else
                switch opts.synthesis.lyapunovshape
                    case 1
                        sGammaw(prob.specs(i).in,prob.specs(i).in) = ssqgamma(j)*eye(length(prob.specs(i).in));
                        sGammaz(prob.specs(i).out,prob.specs(i).out) = ssqgamma(j)*eye(length(prob.specs(i).out));
                    case 2
                        sGammaw(prob.specs(i).in,prob.specs(i).in) = ssqgamma(j)*eye(length(prob.specs(i).in));
                        Gammazcst(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out));
                    case 3
                        Gammawcst(prob.specs(i).in,prob.specs(i).in) = eye(length(prob.specs(i).in));
                        sGammaz(prob.specs(i).out,prob.specs(i).out) = ssqgamma(j)*eye(length(prob.specs(i).out));
                    otherwise
                        prob.watchdog.error('The Lyapunov shaping format can only be of type 1 (-gamma,-gamma), type 2 (-gamma^2,-1) or type 3 (-1,-gamma^2).'); 
                end   
        	end
        end
        
        % transform to LMI variables
        Gammaw = lmivar(3,sGammaw);
        Gammaz = lmivar(3,sGammaz);
        D = lmivar(3,sD);
        
        % constant performance feedthrough matrices
        Dzwcst = prob.Dzwbar();
        Dzucst = prob.Dzubar();
        Dywcst = prob.Dywbar(); 
        Dcst = zeros(size(Dzwcst)); 
        for i=1:length(prob.specs)
            Dcst(prob.specs(i).out,prob.specs(i).in) = Dzwcst(prob.specs(i).out,prob.specs(i).in);
        end
        
    % formulate LMIs
        % STABILITY 
        STAB = newlmi(); 
        lmiterm([-STAB,1,1,X],1,1); 
        lmiterm([-STAB,2,2,Y],1,1); 
        lmiterm([-STAB,2,1,0],prob.XI());
        lmiterm([-STAB,1,1,0],opts.synthesis.zerotol);
        lmiterm([-STAB,2,2,0],opts.synthesis.zerotol);

        % H-INFINITY PERFORMANCE
        PERF = newlmi();
        lmiterm([PERF,1,1,Y],prob.Ahat(),1,'s');
        lmiterm([PERF,1,1,Nc],prob.Buhat(),1,'s');
        lmiterm([PERF,2,1,Kc],1,1);
        lmiterm([PERF,2,1,-Dc],prob.Cytilde()',prob.Buhat()');
        lmiterm([PERF,2,1,0],prob.PoTTASPilmi()');
        lmiterm([PERF,3,1,-Dc],prob.Dywtilde()',prob.Buhat()');
        lmiterm([PERF,3,1,0],prob.Bwlmi()');
        lmiterm([PERF,4,1,Y],prob.Czhat(),1);
        lmiterm([PERF,4,1,Nc],prob.Dzuhat(),1);
        lmiterm([PERF,2,2,X],1,prob.Atilde(),'s');
        lmiterm([PERF,2,2,Mc],1,prob.Cytilde(),'s');
        lmiterm([PERF,3,2,X],prob.Bwtilde()',1);
        lmiterm([PERF,3,2,-Mc],prob.Dywtilde()',1);
        lmiterm([PERF,4,2,Dc],prob.Dzuhat(),prob.Cytilde());
        lmiterm([PERF,4,2,0],prob.Czlmi());
        if ~isempty(Gammaw); lmiterm([PERF,3,3,Gammaw],-1,1); end; lmiterm([PERF,3,3,0],-Gammawcst);
        if ~isempty(Gammaz); lmiterm([PERF,4,4,Gammaz],-1,1); end; lmiterm([PERF,4,4,0],-Gammazcst); 
        if ~isempty(D); lmiterm([PERF,4,3,D],1,1); end; lmiterm([PERF,4,3,0],Dcst);
        Dzumask = zeros(prob.nz(),prob.nu());
        Dywmask = zeros(prob.ny(),prob.nw()); 
        for i=1:length(prob.specs)
            Dzu_ = Dzumask;
            Dyw_ = Dywmask;
            Dzu_(prob.specs(i).out,:) = Dzucst(prob.specs(i).out,:);
            Dyw_(:,prob.specs(i).in) = Dywcst(:,prob.specs(i).in);
            lmiterm([PERF,4,3,Dc],Dzu_,Dyw_);
        end
        lmiterm([PERF,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF,3,3,0],-opts.synthesis.zerotol);
        lmiterm([PERF,4,4,0],-opts.synthesis.zerotol);
        
        % D-STABILITY 
        DSTAB = zeros(length(prob.region),1); 
        for i=1:length(DSTAB)
            if ~(isempty(prob.region(i).M) && isempty(prob.region(i).L))
                DSTAB(i) = newlmi(); 
                for k=1:length(prob.region(i).L)
                    for l=k:length(prob.region(i).L)
                        if l>k
                            % first term
                            lmiterm([DSTAB(i),2*k-1,2*l-1,Y],prob.region(i).L(k,l),1);
                            lmiterm([DSTAB(i),2*k,2*l-1,0],prob.region(i).L(k,l)*prob.XI()');
                            lmiterm([DSTAB(i),2*k,2*l,X],prob.region(i).L(k,l),1);

                            % second term
                            lmiterm([DSTAB(i),2*k-1,2*l-1,Y],prob.region(i).M(k,l)*prob.Ahat(),1);
                            lmiterm([DSTAB(i),2*k-1,2*l-1,Nc],prob.region(i).M(k,l)*prob.Buhat(),1);
                            lmiterm([DSTAB(i),2*k-1,2*l-1,-Y],prob.region(i).M(l,k),prob.Ahat()');
                            lmiterm([DSTAB(i),2*k-1,2*l-1,-Nc],prob.region(i).M(l,k),prob.Buhat()');
                            lmiterm([DSTAB(i),2*k-1,2*l,Dc],prob.region(i).M(k,l)*prob.Buhat(),prob.Cytilde());
                            lmiterm([DSTAB(i),2*k-1,2*l,0],prob.region(i).M(k,l)*prob.PoTTASPilmi());
                            lmiterm([DSTAB(i),2*k-1,2*l,-Kc],prob.region(i).M(l,k),1);
                            lmiterm([DSTAB(i),2*k,2*l,X],prob.region(i).M(k,l),prob.Atilde());
                            lmiterm([DSTAB(i),2*k,2*l,Mc],prob.region(i).M(k,l),prob.Cytilde());
                            lmiterm([DSTAB(i),2*k,2*l,-X],prob.region(i).M(l,k)*prob.Atilde()',1);
                            lmiterm([DSTAB(i),2*k,2*l,-Mc],prob.region(i).M(l,k)*prob.Cytilde()',1);
                        else % l == k -> explicitly tell LMILAB the diagonal block terms are symmetric, otherwise it'll whine
                            % first term
                            lmiterm([DSTAB(i),2*k-1,2*l-1,Y],prob.region(i).L(k,l),1);
                            lmiterm([DSTAB(i),2*k,2*l-1,0],prob.region(i).L(k,l)*prob.XI()');
                            lmiterm([DSTAB(i),2*k,2*l,X],prob.region(i).L(k,l),1);

                            % second term
                            lmiterm([DSTAB(i),2*k-1,2*l-1,Y],prob.region(i).M(k,l)*prob.Ahat(),1,'s');
                            lmiterm([DSTAB(i),2*k-1,2*l-1,Nc],prob.region(i).M(k,l)*prob.Buhat(),1,'s');
                            lmiterm([DSTAB(i),2*k-1,2*l,Dc],prob.region(i).M(k,l)*prob.Buhat(),prob.Cytilde());
                            lmiterm([DSTAB(i),2*k-1,2*l,0],prob.region(i).M(k,l)*prob.PoTTASPilmi());
                            lmiterm([DSTAB(i),2*k-1,2*l,-Kc],prob.region(i).M(l,k),1);
                            lmiterm([DSTAB(i),2*k,2*l,X],prob.region(i).M(k,l),prob.Atilde(),'s');
                            lmiterm([DSTAB(i),2*k,2*l,Mc],prob.region(i).M(k,l),prob.Cytilde(),'s');
                        end

                        lmiterm([DSTAB(i),2*k-1,2*l-1,0],-opts.synthesis.zerotol);
                        lmiterm([DSTAB(i),2*k,2*l,0],-opts.synthesis.zerotol); 
                    end
                end
            end
        end
                
   % optimize
   problem = getlmis(); 
   c = zeros(decnbr(problem),1);
   if any(isobj)
        c(ssqgamma) = prob.specs(isobj).weight; 
   else
        c = [];
   end
   if isnumeric(opts.lmilab) && isvector(opts.lmilab)
       if isempty(c)
           [~,xopt] = feasp(problem,opts.lmilab);
       else
           [~,xopt] = mincx(problem,c,opts.lmilab);
       end
   else
        warning('Unknown options for LMILAB. Ignoring them. Check the mincx() documentation for details.'); 
        [~,xopt] = mincx(problem,c);
   end
 
   % put into the output structure
   sdpvars.X = dec2mat(problem,xopt,X);
   sdpvars.Y = dec2mat(problem,xopt,Y); 
   sdpvars.Kc = dec2mat(problem,xopt,Kc); 
   sdpvars.Mc = dec2mat(problem,xopt,Mc);
   sdpvars.Nc = dec2mat(problem,xopt,Nc);
   sdpvars.Dc = dec2mat(problem,xopt,Dc); 
   sdpvars.sqgamma = dec2mat(problem,xopt,sqgamma); 
  
end

%% CONTINUOUS TIME - RECONSTRUCTION

function K = ct_controller_reconstruction(sdpvars,prob)

    % do the first backtransformation
    Ka = sdpvars.Kc - sdpvars.X(:,1:prob.n())*prob.A()*sdpvars.Y((prob.no()+1):end,:);
    Ma = sdpvars.Mc;
    Na = sdpvars.Nc; 
    Da = sdpvars.Dc; 
    
    % do the second backtransformation
    [PiTU,S,PoTV] = svd(prob.XI()' - sdpvars.X(:,1:prob.n())*sdpvars.Y((prob.no()+1):end,:));
    PiTU = PiTU*sqrt(S(:,1:prob.n()));
    PoTV = PoTV*sqrt(S(1:prob.n(),:))'; 
    
    PiTXSiBe = -sdpvars.X(:,1:prob.n())*[prob.PIi() -prob.Bu()] + [sdpvars.X(:,(prob.n()+1):end) zeros(prob.n()+prob.ni(),prob.nu())];
    CeTiYPo = [-prob.PIo() ; prob.Cy()]*sdpvars.Y((prob.no()+1):end,:) + [sdpvars.Y(1:prob.no(),:) ; zeros(prob.ny(),prob.n()+prob.no())]; 
    PRE = [PiTU PiTXSiBe ; zeros(prob.nu(),prob.n()+prob.ni()) eye(prob.nu())];
    POST = [[PoTV' ; CeTiYPo] [zeros(prob.n()+prob.no(),prob.ny()) ; eye(prob.ny())]]; 
    
    THETA = PRE\[Ka Ma ; Na Da]/POST; 
    
    % cast state-space controller matrices
    Ak = THETA(1:prob.n(),1:prob.n());
    Bk = THETA(1:prob.n(),(prob.n()+1):end);
    Ck = THETA((prob.n()+1):end,1:prob.n()); 
    Dk = THETA((prob.n()+1):end,(prob.n()+1):end); 
    
    % return the controller
    K = ss(Ak,Bk,Ck,Dk);

end
