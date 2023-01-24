%% Source file of design example - plain sshinfcd
%
%  See Section V of Jacobs & Swevers (2023). "H control design with
%  D-unstable weighting filters and a D-stability constraint: Solution and
%  applications." Submitted for the 7th IEEE Conference on Control
%  Technology and Applications (CCTA). 
%
%  For this specific problem, we use CVX (http://cvxr.com/cvx/) and 
%  MOSEK (https://www.mosek.com/). 

% This file is part of sshinfcd.
% Copyright (c) 2023, Laurens Jacobs, MECO Research Team @ KU Leuven. 
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
pmodel = canon(fit,'modal');
npmodel = G;

%% construct the generalized plant
%  from (reference, input disturbance, control signal) -> (error, output, control signal, error)
GP = minreal([1 -pmodel -pmodel;
              0  pmodel  pmodel; 
              0       0       1;
              1 -pmodel -pmodel]); 
          
%% define the weighting filters and the specs

% individual weights
WS = ss(0,4,5.230751768227005,0);                                           % low-frequency roll-off 
MS = 0.316227766016838;                                                     % peak gain
WU = ss([-18856.0601676690, -15233.8254177579, -5859.69653074140;
          15233.8254177579, -2274.31914203070, -3859.23561308874;
         -5859.69653074140,  3859.23561308875, -8033.09736532489], ...
        [ 1729.46680178491;
         -395.583915265197;
          285.970073684068], ...
        [-1729.46680178491, -395.583915265195, -285.970073684068], ...
        100);                                                               % high-frequency roll-off
 
% concatenate
Wi = ss(eye(2));
Wo = [WS  0   0; 
      MS  0   0; 
       0  W   0; 
       0  0  WU];

% define the pole region constraint
zeta = 0.2;                                                                 % minimal damping ratio
alpha = zeros(2); 
beta = [sqrt(1-zeta^2)               -zeta;
        zeta                sqrt(1-zeta^2)]; 

% define the (mixed-sensitivity) specs
specs = struct(); 
specs(1).in = 1; specs(1).out = 1; specs(1).weight = 0; 
specs(2).in = 1; specs(2).out = 2; specs(2).weight = 0;
specs(3).in = 1; specs(3).out = 3; specs(3).weight = 0;
specs(4).in = 1; specs(4).out = 4; specs(4).weight = 1;

% define the pole region constraints
reg(1).L = alpha; 
reg(1).M = beta;

%% create the solver and set options

solver = sshinfcd.sshinfcd(); 
opts = solver.options(); 
opts.synthesis.lyapunovshape = 3; 
opts.synthesis.gammasolver = 'cvx';
opts.cvx.solver = 'mosek'; 

%% solve the controller design problem without damping

solver = solver.setproblem(GP,Wi,Wo,specs);
[K1,gamma1] = solver.solve(opts); 

%% solve the controller design problem with damping

solver = solver.setproblem(GP,Wi,Wo,specs,reg);
[K2,gamma2] = solver.solve(opts); 

%% compare results

CL1 = lft(GP,K1); 
CL2 = lft(GP,K2); 
figure; bodemag(CL1(1,1),CL2(1,1)); title('Sensitivity'); hold on; bodemag(1/WS,'k--'); bodemag(1/ss(MS),'k--'); legend('without D-stability constraint','with D-stability constraint','',''); 
figure; bodemag(CL1(2,1),CL2(2,1)); title('Complementary sensitivity'); hold on; bodemag(1/W,'k--'); legend('without D-stability constraint','with D-stability constraint',''); 
figure; bodemag(CL1(3,1),CL2(3,1)); title('Control sensitivity'); hold on; bodemag(gamma1/WU,'k--'); bodemag(gamma2/WU,'k--'); legend('without D-stability constraint','with D-stability constraint',''); 

p1 = pole(CL1);
p2 = pole(CL2);
figure; plot(real(p1),imag(p1),'x'); hold on; plot(real(p2),imag(p2),'x'); title('Pole-zero map'); 
xlabel('Real (rad/s)'); ylabel('Im (rad/s)'); xlim([-600 0]); ylim([-1000 1000]); 
plot([-204.12 0 -204.12],[-1000 0 1000],'k--'); 