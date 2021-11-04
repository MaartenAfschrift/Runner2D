% Default example 2D predictive simulations
%-------------------------------------------

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1;        % average speed
S.N         = 140;       % number of mesh intervals

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.93;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'Example_2DSim';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Gait18_subject1'; 

% select folder with polynomials
S.PolyFolder = 'Gait18_subject1';

% initial guess based on simulations without exoskeletons
S.IDig          = 2;        % (1) == quasi random, (2) == IK file, (3) from solution
% S.ResultsF_ig   = 'Default';
% S.savename_ig   = '_a_SimMTP_85_q';

% Select model
S.ModelName = 'Gait18'; % other option is 'Rajagopal'

% mass of the subject
S.mass = 64;    % mass in kg

% IK folder
% S.IKFolder = 's1_Poggensee';

% linear solver
S.linear_solver = 'ma57';

% Simulation without exoskeelton
S.ExternalFunc = 'SimExo_mtp_2D';
% S.ExternalFunc  = 'PredSim_2D_Pog_s1_mtp.dll';
S.savename      = 'DefaultSim_2D_v5';

% predictive simulation
f_PredSim_2D_euler(S);

% CUrrent state: I get NaN when using metabolic energy in the objective
% function. Origin of these NaNs is unclear !