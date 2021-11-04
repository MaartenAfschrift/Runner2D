% Default example 2D predictive simulations
%-------------------------------------------

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.3;        % average speed
S.N         = 50;       % number of mesh intervals

% path repository
S.pathRepo = 'C:\Users\u0088756\Documents\Software\PredSim_2D';

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.93;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'Sensitivity_WeightsObj2';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Gait18_shorterHamstring';

% select folder with polynomials
S.PolyFolder ='Gait18_shorterHamstring';

% initial guess based on simulations without exoskeletons
S.IDig          = 2;        % (1) == quasi random, (2) == IK fisdle, (3) from solution
S.ResultsF_ig   = 'Example_2DSim';
S.savename_ig   = 'DefaultSim_2D_trapz3_q';

% Select model
S.ModelName = 'Gait18'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg

% linear solver
S.linear_solver = 'mumps';

% Simulation without exoskeelton
% S.ExternalFunc = 'SimExo_mtp_2D';
S.ExternalFunc  = 'PredSim_2D_default_mtp_v3.dll';
S.ExternalFunc_pp  = 'PredSim_2D_default_mtp_v3_pp.dll';
S.savename      = S.CasadiFunc_Folders; 

% set the weights in the objective function
S.W.A         = 2;    % weight muscle activations
S.W.back      = 1;    % weight back excitations
S.W.acc       = 2;    % weight joint accelerations
S.W.exp_A     = 3;    % power muscle activations
S.W.u         = 0.001;    % power muscle activations
S.W.E         = 0;    % metabolic energy

% tolerance ipopt
S.tol_ipopt   = 6;    % tolerance ipopt

% open a cluster
myCluster = parcluster('Maarten_LocalProfile1'); % you can create your own cluster profile...

% first test some sensitivities from current default parameters
AVect = [0 0.01: 0.1:1 1:0.5:3 3:2:10 50 100];
AccVect = [0 0.01: 0.1:1 1:0.5:3 3:2:10 50 100];
EVect = [0 0.001:0.01:0.05 0.1:0.2:2 50 100];

ct = 1;
for i=1:length(AVect)
    S.W.A       = AVect(i);        % weight muscle activations    
    S.W.back      = 1;    % weight back excitations
    S.W.acc       = 2;    % weight joint accelerations
    S.W.exp_A     = 3;    % power muscle activations
    S.W.u         = 0.001;    % power muscle activations
    S.W.E         = 0;    % metabolic energy
    S.savename  = ['SensObj_2D_wA_' num2str(i)];
    batch(myCluster,'f_PredSim_2D_trapezoidal',0,{S});
    ct = ct+1;
end

for i=1:length(AccVect)
    S.W.A         = 2;    % weight muscle activations
    S.W.back      = 1;    % weight back excitations
    S.W.acc       = AccVect(i);    % weight joint accelerations
    S.W.exp_A     = 3;    % power muscle activations
    S.W.u         = 0.001;    % power muscle activations
    S.W.E         = 0;    % metabolic energy
    S.savename  = ['SensObj_2D_wAcc_' num2str(i)];
    batch(myCluster,'f_PredSim_2D_trapezoidal',0,{S});
    ct = ct+1;
end

Nrand = 5000;
Wrand = randn(Nrand,3)*100;
for i=1:Nrand
    S.W.A         = Wrand(i,1);    % weight muscle activations
    S.W.back      = Wrand(i,2);    % weight back excitations
    S.W.acc       = Wrand(i,3);    % weight joint accelerations
    S.W.exp_A     = 3;    % power muscle activations
    S.W.u         = 0.001;    % power muscle activations
    S.W.E         = 0;    % metabolic energy
    S.savename  = ['SensObj_2D_wRand_' num2str(i)];
    batch(myCluster,'f_PredSim_2D_trapezoidal',0,{S});
    ct = ct+1;
end