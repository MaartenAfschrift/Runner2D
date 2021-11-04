% Default example 2D predictive simulations
%-------------------------------------------

clear all; close all; clc;
% settings for optimization
% S.v_tgt     = 1.3;        % average speed
S.N         = 50;       % number of mesh intervals

% path repository
S.pathRepo = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D';

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.93;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'Example_OptSpeed';

% selection folder with Casadi Functions
S.CasadiFunc_Folders = 'Gait18_subject1'; 

% select folder with polynomials
S.PolyFolder = 'Gait18_subject1';

% initial guess based on simulations without exoskeletons
S.IDig          = 2;        % (1) == quasi random, (2) == IK file, (3) from solution
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
S.ExternalFunc  = 'FastRunner.dll';
S.savename      = 'testFast';

S.k_mtp = 100;
% S.d_mtp = 1;
% predictive simulation
[Results] = f_PredSim_2D_trapezoidal_OptSpeed(S);

%% visualise the solution
% 
% ModelPath = fullfile(S.pathRepo,'OpenSimModel','Gait18_Antoine2.osim'); 
% motFile = fullfile(S.pathRepo,'Results',S.ResultsFolder,[S.savename '_q.mot']);
% Visualise2DModel(ModelPath,motFile,2);


%% Plot some results
% figure(); 
% plot([Results.lMtilde(:,1:9);Results.lMtilde(:,10:18)]); 
% legend(Results.muscleNames);
% 
% figure(); 
% plot([Results.a(:,1:9);Results.a(:,10:18)]); 
% legend(Results.muscleNames);

disp(['Running speed is ' num2str(Results.speed*3.6)]);



