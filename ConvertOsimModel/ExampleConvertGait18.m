% Example script convert gait18 model

% path information
MainPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D';
% settings:
% Folder to save the polynomials
S.PolyFolder = 'Gait18_shorterHamstring';
% Modelpath
S.ModelPath = fullfile(MainPath,'OpenSimModel','Gait18_shorterHamstring.osim'); 
% Folder with CasadiFunctions
S.CasadiFunc_Folders = 'Gait18_shorterHamstring'; 
% model selection options: Rajagopal, Gait92
S.ModelName = 'Gait18';      

% specific settings for exporting casadi functions
SettingsCasFunc.kTendon_CalfM = 20;
SettingsCasFunc.kMTP = 1.5/(pi/180)/5;
SettingsCasFunc.dMTP = 0.5;
%% Fit polynmial functions

Bool_RunMA = true; % Boolean to select if we have to run the muscle analysis
FitPolynomials(MainPath,S.ModelName,S.ModelPath,...
    S.PolyFolder,Bool_RunMA);


%% Create only the casadi functions
% create casadi functions for equations in optimiztion problem
CreateCasadiFunctions(MainPath, S.ModelName, S.ModelPath, S.CasadiFunc_Folders,...
    S.PolyFolder,SettingsCasFunc);

