%% Source file of design example - with LCToolbox + sshinfcd
%
%  See Section V of Jacobs & Swevers (2023). "H control design with
%  D-unstable weighting filters and a D-stability constraint: Solution and
%  applications." Submitted for the 7th IEEE Conference on Control
%  Technology and Applications (CCTA). 
%
%  To run this file, you need LCToolbox, which is available online: 
%  https://github.com/meco-group/lc_toolbox
%
%  For this specific problem, we also use CVX (http://cvxr.com/cvx/) and 
%  MOSEK (https://www.mosek.com/). 

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

clear; close all; 

% load design models
load('model.mat');
pmodel = fromstd(canon(fit,'modal'));
pmodel.name = 'parametric model'; 
npmodel = fromstd(G);
npmodel.name = 'nonparametric model'; 

%% controller design with damping

% weights
MS = Weight.DC(10);
WS = 3.33*Weight.LF(1,1);
WU = Weight.HF(500,3,-40);
 
params.beta = asin(0.2); 
reg = Region('conic',params); 

P = IOSystem(1,1);
P.add(pmodel); 
P.add(npmodel); 
K = IOSystem(1,1);

r = Signal(1);      % reference
d = Signal(1);      % disturbance
u = K.out;          % control signal
y = P.out;          % response
e = r-y;            % tracking error

conn = [K.in == e; 
        P.in == u+d];
CL = IOSystem(P,K,conn);

S = Channel(e/r,'S - Sensitivity');
Ud = Channel(u/d,'Ud - Input sensitivity');
Ur = Channel(u/r,'Ur - Input sensitivity');
T = Channel(y/r,'T - Complementary sensitivity');
D = Channel(y/d,'D - Disturbance sensitivity'); 

obj = [WU*Ur];     
cstr = [WS*S <= 1;
        MS*S <= 1;
        W*T <= 1; 
        reg];       
opts.FullOrderSolver = 'sshinfcd';
opts.synthesis.lyapunovshape = 3;
opts.synthesis.gammasolver = 'cvx';
opts.cvx.solver = 'mosek'; 
opts.controller_name = 'with D-stability constraint';
[~,~,info] = CL.solve(obj,cstr,K,opts); 

%% controller design without damping

% weights (= same as above)
MS = Weight.DC(10);
WS = 3.33*Weight.LF(1,1);
WU = Weight.HF(500,3,-30);

obj = [WU*Ur];     
cstr = [MS*S <= 1;
        W*T <= 1; 
        WS*S <= 1];
opts.FullOrderSolver = 'sshinfcd';
opts.synthesis.lyapunovshape = 3;
opts.synthesis.gammasolver = 'cvx';
opts.cvx.solver = 'mosek'; 
opts.controller_name = 'without D-stability constraint';
[~,~,info2] = CL.solve(obj,cstr,K,opts); 

%% compare results

showall(info,info2);
figure; bode(K); 

CL = CL(S); 
CL1 = CL.content(1);
CL2 = CL.content(2);
p1 = pole(CL1);
p2 = pole(CL2);
figure; plot(real(p1),imag(p1),'x'); hold on; plot(real(p2),imag(p2),'x'); title('Pole-zero map'); 
xlabel('Real (rad/s)'); ylabel('Im (rad/s)'); xlim([-600 0]); ylim([-1000 1000]); 
plot([-204.12 0 -204.12],[-1000 0 1000],'k--'); 