%% Casadi Functions for energy models
%-------------------------------------

clear all; close all; clc;

import casadi.*

OutPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D\CasadiFunctions\EnergyModels';
if ~isfolder(OutPath)
    mkdir(OutPath);
end

%% Casadi variables
NMuscle         = 18;
act_SX          = SX.sym('act_SX',NMuscle,1); % Muscle activations
exc_SX          = SX.sym('exc_SX',NMuscle,1); % Muscle excitations
lMtilde_SX      = SX.sym('lMtilde_SX',NMuscle,1); % N muscle fiber lengths
vMtilde_SX      = SX.sym('vMtilde_SX',NMuscle,1); % N muscle fiber vel
vM_SX           = SX.sym('vM_SX',NMuscle,1); % Muscle fiber velocities
Fce_SX          = SX.sym('FT_SX',NMuscle,1); % Contractile element forces
Fpass_SX        = SX.sym('FT_SX',NMuscle,1); % Passive element forces
Fiso_SX         = SX.sym('Fiso_SX',NMuscle,1); % N forces (F-L curve)
musclemass_SX   = SX.sym('musclemass_SX',NMuscle,1); % Muscle mass
vcemax_SX       = SX.sym('vcemax_SX',NMuscle,1); % Max contraction vel
pctst_SX        = SX.sym('pctst_SX',NMuscle,1); % Slow twitch ratio
Fmax_SX         = SX.sym('Fmax_SX',NMuscle,1); % Max iso forces
modelmass_SX    = SX.sym('modelmass_SX',1); % Model mass
b_SX            = SX.sym('b_SX',1); % Parameter determining tanh smoothness



%% Energy models

% Bhargava et al. (2004)
[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2004all(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,pctst_SX,Fiso_SX,...
    Fmax_SX,modelmass_SX,b_SX);
fgetMetabolicEnergySmooth2004all = Function('fgetMetabolicEnergySmooth2004all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,...
    pctst_SX,Fiso_SX,Fmax_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_model_sm_SX},...
    {'e','a','lMtilde','vM','Fce','Fpass','musclemass',...
    'pctst','Fiso','Fmax','modelmass','b'},...
    {'energy_total','Adot','Mdot','Sdot','Wdot','energy_model'});

% Umberger 2003
[energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2003all(exc_SX,act_SX,...
    lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,pctst_SX,...
    vcemax_SX,Fiso_SX,modelmass_SX,b_SX);                                
fgetMetabolicEnergySmooth2003all = Function('fgetMetabolicEnergySmooth2003all',...
    {exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,...
    pctst_SX,vcemax_SX,Fiso_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,energy_model_sm_SX},...
    {'e','a','lMtilde','vMtilde','vM','Fce','musclemass','pctst','vcemax','Fiso','modelmass','b'},...
    {'energy_total','energy_am','energy_sl','energy_mech','energy_model'});

% Umberger 2010
[energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2010all(exc_SX,act_SX,...
    lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,pctst_SX,vcemax_SX,...
    Fiso_SX,modelmass_SX,b_SX);                                
fgetMetabolicEnergySmooth2010all = ...
    Function('fgetMetabolicEnergySmooth2010all',...
    {exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,...
    pctst_SX,vcemax_SX,Fiso_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,energy_model_sm_SX},...
    {'e','a','lMtilde','vMtilde','vM','Fce','musclemass','pctst','vcemax','Fiso','modelmass','b'},...
    {'energy_total','energy_am','energy_sl','energy_mech','energy_model'});

% Uchida 2016
[energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2016all(exc_SX,act_SX,...
    lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,pctst_SX,vcemax_SX,...
    Fiso_SX,modelmass_SX,b_SX);                                
fgetMetabolicEnergySmooth2016all = ...
    Function('fgetMetabolicEnergySmooth2016all',...
    {exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,...
    pctst_SX,vcemax_SX,Fiso_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,energy_model_sm_SX},...
    {'e','a','lMtilde','vMtilde','vM','Fce','musclemass','pctst','vcemax','Fiso','modelmass','b'},...
    {'energy_total','energy_am','energy_sl','energy_mech','energy_model'});  

% Umberger et al. (2003)                       
[energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2010all_hl(exc_SX,...
    act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,pctst_SX,...
    vcemax_SX,Fiso_SX,modelmass_SX,b_SX);                                
fgetMetabolicEnergySmooth2010all_hl = ...
    Function('fgetMetabolicEnergySmooth2010all_hl',...
    {exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,...
    pctst_SX,vcemax_SX,Fiso_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,energy_model_sm_SX},...
    {'e','a','lMtilde','vMtilde','vM','Fce','musclemass','pctst','vcemax','Fiso','modelmass','b'},...
    {'energy_total','energy_am','energy_sl','energy_mech','energy_model'});

% Umberger et al. (2003)
[energy_total_sm_SX,energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2010all_neg(exc_SX,...
    act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,pctst_SX,...
    vcemax_SX,Fiso_SX,modelmass_SX,b_SX);                                
fgetMetabolicEnergySmooth2010all_neg = ...
    Function('fgetMetabolicEnergySmooth2010all_neg',...
    {exc_SX,act_SX,lMtilde_SX,vMtilde_SX,vM_SX,Fce_SX,musclemass_SX,...
    pctst_SX,vcemax_SX,Fiso_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    energy_am_sm_SX,energy_sl_sm_SX,energy_mech_sm_SX,energy_model_sm_SX},...
    {'e','a','lMtilde','vMtilde','vM','Fce','musclemass','pctst','vcemax','Fiso','modelmass','b'},...
    {'energy_total','energy_am','energy_sl','energy_mech','energy_model'});


% energy equations Alexander



% energy equations Minetti



% energy equations Margaria:
[metabE] = getMetabolicEnergy_MargariaSmooth(Fce_SX,vM_SX,b_SX);
fgetMetabolicEnergy_MargariaSmooth = Function('fgetMetabolicEnergy_MargariaSmooth',...
    {Fce_SX,vM_SX,b_SX},{metabE},{'Fce','vM','b'},{'Edot'});




%% save the casadi functions

% export the casadi functions
fgetMetabolicEnergySmooth2004all.save(fullfile(OutPath,'fgetMetabolicEnergySmooth2004all'));
fgetMetabolicEnergySmooth2003all.save(fullfile(OutPath,'fgetMetabolicEnergySmooth2003all'));
fgetMetabolicEnergySmooth2010all.save(fullfile(OutPath,'fgetMetabolicEnergySmooth2010all'));
fgetMetabolicEnergySmooth2016all.save(fullfile(OutPath,'fgetMetabolicEnergySmooth2016all'));
fgetMetabolicEnergySmooth2010all_hl.save(fullfile(OutPath,'fgetMetabolicEnergySmooth2010all_hl'));
fgetMetabolicEnergySmooth2010all_neg.save(fullfile(OutPath,'fgetMetabolicEnergySmooth2010all_neg'));
fgetMetabolicEnergy_MargariaSmooth.save(fullfile(OutPath,'fgetMetabolicEnergy_MargariaSmooth'));