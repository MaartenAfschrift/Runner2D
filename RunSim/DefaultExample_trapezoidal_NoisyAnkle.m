% Default example 2D predictive simulations
%-------------------------------------------

clear all; close all; clc;
% settings for optimization
S.v_tgt     = 1.3;        % average speed
S.N         = 50;       % number of mesh intervals

% path repository
S.pathRepo = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D';

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.93;   % subject 1 poggensee

% output folder
S.ResultsFolder     = 'Example_2DSim_NoisyAnkle';

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

% predictive simulation
[Results] = f_PredSim_2D_trapezoidal_NoisyAnkle(S);

%% visualise the solution

% ModelPath = fullfile(S.pathRepo,'OpenSimModel','Gait18_Visual.osim'); 
% motFile = fullfile(S.pathRepo,'Results',S.ResultsFolder,[S.savename '_q.mot']);
% Visualise2DModel(ModelPath,motFile,2);

% figure(); 
% plot([Results.lMtilde(:,1:9);Results.lMtilde(:,10:18)]); 
% legend(Results.muscleNames);
% 
% figure(); 
% plot([Results.a(:,1:9);Results.a(:,10:18)]); 
% legend(Results.muscleNames);
% 
% figure();
% for i=1:3
%    subplot(1,3,i)
%    MA = squeeze(Results.MA(:,:,i));   
%    plot([MA(:,1:9); MA(:,10:18)]);    
%    legend(Results.muscleNames);
% 
% end
