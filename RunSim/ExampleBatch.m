% Default example 2D predictive simulations
%-------------------------------------------

clear all; close all; clc;

% settings for optimization
S.N         = 50;       % number of mesh intervals

% main path repository
S.pathRepo = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D';

% output folder
S.ResultsFolder     = 'Example_Batch';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Gait18_subject1'; 

% select folder with polynomials
S.PolyFolder = 'Gait18_subject1';

% initial guess based on simulations without exoskeletons
S.IDig          = 3;        % (1) == quasi random, (2) == IK file, (3) from solution
S.ResultsF_ig   = 'Example_2DSim';
S.savename_ig   = 'DefaultSim_2D_trapz3_q';

% Select model
S.ModelName = 'Gait18'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg

% linear solver
S.linear_solver = 'mumps';

% Simulation without exoskeelton
S.ExternalFunc = 'SimExo_mtp_2D';

% predictive simulation
Speeds = 0.5:0.1:1.3;

myCluster = parcluster('local'); % you can create your own cluster profile...

for i=1:length(Speeds)
    S.v_tgt     = Speeds(i);        % average speed
    S.savename  = ['DefaultSim_2D_speed_' num2str(round(Speeds(i)*10))];
    Jobs(i) = batch(myCluster,'f_PredSim_2D_trapezoidal',0,{S});
end
