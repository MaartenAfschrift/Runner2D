function [] = CreateCasadiFunctions(MainPath, ModelName, ModelPath, CasadiFunc_Folders,...
    PolyFolder,Settings,varargin)

%CreateCasadiFunctions Function to create all casadifunctions needed in the
%optimization

import casadi.*

boolOverWrite = false;
if ~isempty(varargin)
    boolOverWrite = varargin{1};
end

if isfolder(fullfile(MainPath,'CasADiFunctions',CasadiFunc_Folders)) && ~boolOverWrite
    error(['Never changes the casadi functions in an existing folder,',...
        'these functions are important to analyse optimization results',...
        '(from the optimal states and controls. Please change the name of',...
        ' the folder with casadifunctions']);
end


%% update settings structure with defaults
if ~exist('Settings','var')
    Settings.kTendon_CalfM = 20;        % Atendon of calf muscles
    Settings.kMTP = 1.5/(pi/180)/5;     % stiffness mtp joint
    Settings.dMTP = 0.5;                % damping mtp joint
else
    if ~isfield(Settings,'kTendon_CalfM')
        Settings.kTendon_CalfM = 20;        % Atendon of calf muscles
    end
    if ~isfield(Settings,'kMTP')
        Settings.kMTP = 1.5/(pi/180)/5;
    end
    if ~isfield(Settings,'dMTP')
        Settings.dMTP = 0.5;
    end
end

%% specific settings depending on the model used
if strcmp(ModelName,'Gait18')
    % Indices of the elements in the external functions
    jointi.pelvis.tilt  = 1;
    jointi.pelvis.tx    = 2;
    jointi.pelvis.ty    = 3;
    jointi.hip.l        = 4;
    jointi.hip.r        = 5;
    jointi.knee.l       = 6;
    jointi.knee.r       = 7;
    jointi.ankle.l      = 8;
    jointi.ankle.r      = 9;
    jointi.mtp.l        = 10;
    jointi.mtp.r        = 11;
    jointi.trunk.ext    = 12;
    % Vectors of indices for later use
    jointi.all          = jointi.pelvis.tilt:jointi.trunk.ext; % all
    jointi.gr_pelvis    = jointi.pelvis.tilt:jointi.pelvis.ty; % ground-pelvis
    % Number of degrees of freedom for later use
    nq.all              = length(jointi.all); % all
    nq.abs              = length(jointi.gr_pelvis); % ground-pelvis
    nq.leg              = 3; % #joints needed for polynomials
    nq.trunk            = 1; % trunk
    % External function: F1 (post-processing purpose only)
    % Ground reaction forces (GRFs)
    GRFi.r              = 13:14;
    GRFi.l              = 15:16;
    GRFi.all            = [GRFi.r,GRFi.l];
    
    % Muscle indices for later use
    muscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
        'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};
    musi    = MuscleIndices_2D(muscleNames);
    NMuscle = length(muscleNames)*2;
end
    
    % Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation
% angles; Row 5: maximal contraction velocities
MTparameters =  ReadMuscleParameters(ModelPath,muscleNames);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];

% Indices of the muscles actuating the different joints for later use
MinfoSpanning = load(fullfile(MainPath,'Polynomials',PolyFolder,'muscle_spanning_joint_INFO.mat'));
[Indmusi,mai] = MomentArmIndices_2D(muscleNames,MinfoSpanning.muscle_spanning_joint_INFO);

% adapt the tendon stiffness
TendonInfo.Atendon = ones(NMuscle,1)*35;
TendonInfo.Atendon([7:8 16:17]) = Settings.kTendon_CalfM;
TendonInfo.shift = getShift(TendonInfo.Atendon);

%% Metabolic energy model parameters

% We load some variables for the polynomial approximations
MInfo = load(fullfile(MainPath,'Polynomials',PolyFolder,'MuscleInfo.mat'));
musi_pol = musi;
NMuscle_pol = NMuscle/2;


%% polynomial equations
muscle_spanning_info_m = MinfoSpanning.muscle_spanning_joint_INFO(musi_pol,:);
MuscleInfo_m.muscle    = MInfo.MuscleInfo.muscle(musi_pol);                  
q     = SX.sym('q',1,nq.leg);
qd    = SX.sym('qd',1,nq.leg);
lMT   = SX(NMuscle_pol,1);
vMT   = SX(NMuscle_pol,1);
dM    = SX(NMuscle_pol,nq.leg);
for i=1:NMuscle_pol      
    index_dof_crossing  = find(muscle_spanning_info_m(i,:)==1);
    order               = MuscleInfo_m.muscle(i).order;
    [mat,diff_mat_q]    = n_art_mat_3_cas_SX(q(1,index_dof_crossing),order);
    lMT(i,1)            = mat*MuscleInfo_m.muscle(i).coeff;
    vMT(i,1)            = 0;
    dM(i,1:nq.leg)      = 0;
    nr_dof_crossing     = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = ...
            (-(diff_mat_q(:,dof_nr)))'*MuscleInfo_m.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*...
            qd(1,index_dof_crossing(dof_nr)));
    end 
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{q,qd},{lMT,vMT,dM},...
    {'q','qd'},{'lMT','vMT','dM'});

%% Muscle contraction dynamics

% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde',NMuscle); % Normalized tendon forces
a           = SX.sym('a',NMuscle); % Muscle activations
dFTtilde    = SX.sym('dFTtilde',NMuscle); % Time derivative tendon forces
lMT         = SX.sym('lMT',NMuscle); % Muscle-tendon lengths
vMT         = SX.sym('vMT',NMuscle); % Muscle-tendon velocities
Hilldiff    = SX(NMuscle,1); % Hill-equilibrium
FT          = SX(NMuscle,1); % Tendon forces
Fce         = SX(NMuscle,1); % Contractile element forces
Fpass       = SX(NMuscle,1); % Passive element forces
Fiso        = SX(NMuscle,1); % Normalized forces from force-length curve
vMmax       = SX(NMuscle,1); % Maximum contraction velocities
massM       = SX(NMuscle,1); % Muscle mass
% Parameters of force-length-velocity curves
load(fullfile(MainPath,'MuscleModel','Fvparam.mat'),'Fvparam');
load(fullfile(MainPath,'MuscleModel','Fpparam.mat'),'Fpparam');
load(fullfile(MainPath,'MuscleModel','Faparam.mat'),'Faparam');

if ~exist('TendonInfo','var')
    TendonInfo.Atendon = ones(NMuscle,1)*35;
    TendonInfo.shift = zeros(NMuscle,1);
end

for m = 1:NMuscle    
    [Hilldiff(m),FT(m),Fce(m),Fpass(m),Fiso(m),vMmax(m)] = ...
        ForceEquilibrium_FtildeState(a(m),FTtilde(m),dFTtilde(m),...
        lMT(m),vMT(m),MTparameters_m(:,m),Fvparam,Fpparam,Faparam,...
        TendonInfo.Atendon(m),TendonInfo.shift(m));    
end
f_forceEquilibrium_FtildeState = ...
    Function('f_forceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,...
        lMT,vMT},{Hilldiff,FT,Fce,Fpass,Fiso,vMmax,massM});
    
% Function to get (normalized) muscle fiber lengths
lM      = SX(NMuscle,1);
lMtilde = SX(NMuscle,1);
for m = 1:NMuscle
    [lM(m),lMtilde(m)] = FiberLength_TendonForce(FTtilde(m),...
        MTparameters_m(:,m),lMT(m),TendonInfo.Atendon(m),TendonInfo.shift(m));
end
f_FiberLength_TendonForce = Function('f_FiberLength_Ftilde',...
    {FTtilde,lMT},{lM,lMtilde});

% Function to get (normalized) muscle fiber velocities
vM      = SX(NMuscle,1);
vMtilde = SX(NMuscle,1);
for m = 1:NMuscle
    [vM(m),vMtilde(m)] = FiberVelocity_TendonForce(FTtilde(m),...
        dFTtilde(m),MTparameters_m(:,m),lMT(m),vMT(m),TendonInfo.Atendon(m),TendonInfo.shift(m));
end
f_FiberVelocity_TendonForce = Function('f_FiberVelocity_Ftilde',...
    {FTtilde,dFTtilde,lMT,vMT},{vM,vMtilde});

%% Back activation dynamics
e_b = SX.sym('e_b',nq.trunk); % back excitations
a_b = SX.sym('a_b',nq.trunk); % back activations
dadt = BackActivationDynamics(e_b,a_b);
f_BackActivationDynamics = ...
    Function('f_BackActivationDynamics',{e_b,a_b},{dadt},...
    {'e_b','a_b'},{'dadt'});

%% Default activation dynamics

e_m = SX.sym('e_b',NMuscle,1); % back excitations
a_m = SX.sym('a_m',NMuscle,1); % back activations
dadt_m = BackActivationDynamics(e_m,a_m);
f_ActivationDynamics = ...
    Function('f_ActivationDynamics',{e_m,a_m},{dadt_m},...
    {'e','a'},{'dadt'});

%% small functions

% Normalized sum of values to a certain power
etempNMuscle = SX.sym('etempNMuscle',NMuscle);
exp          = SX.sym('exp',1);
JtempNMuscle = 0;
for i=1:length(etempNMuscle)
    JtempNMuscle = JtempNMuscle + etempNMuscle(i).^exp;
end
f_sumsqr_exp=Function('f_sumsqr_exp',{etempNMuscle,exp},{JtempNMuscle});

%% Sum of products  
% Function for 3 elements 
ma_temp3 = SX.sym('ma_temp3',3);
ft_temp3 = SX.sym('ft_temp3',3);
J_sptemp3 = 0;
for i=1:length(ma_temp3)
    J_sptemp3 = J_sptemp3 + ma_temp3(i,1)*ft_temp3(i,1);    
end
f_T3 = Function('f_T3',{ma_temp3,ft_temp3},{J_sptemp3});
% Function for 4 elements 
ma_temp4 = SX.sym('ma_temp4',4);
ft_temp4 = SX.sym('ft_temp4',4);
J_sptemp4 = 0;
for i=1:length(ma_temp4)
    J_sptemp4 = J_sptemp4 + ma_temp4(i,1)*ft_temp4(i,1);    
end
f_T4 = Function('f_T4',{ma_temp4,ft_temp4},{J_sptemp4});
% Function for 5 elements 
ma_temp5 = SX.sym('ma_temp5',5);
ft_temp5 = SX.sym('ft_temp5',5);
J_sptemp5 = 0;
for i=1:length(ma_temp5)
    J_sptemp5 = J_sptemp5 + ma_temp5(i,1)*ft_temp5(i,1);    
end
f_T5 = Function('f_T5',{ma_temp5,ft_temp5},{J_sptemp5});

%% activation dynamics
e_m = SX.sym('e_b',NMuscle,1); % back excitations
a_m = SX.sym('a_m',NMuscle,1); % back activations
dadt_m = DefaultActivationDynamics(e_m,a_m);
f_ActivationDynamics = ...
    Function('f_ActivationDynamics',{e_m,a_m},{dadt_m},...
    {'e','a'},{'dadt'});

%% save the casadifunctions

OutPath = fullfile(MainPath,'CasADiFunctions',CasadiFunc_Folders);
if ~isfolder(OutPath)
    mkdir(OutPath);
end
save(fullfile(OutPath,'MTparameters.mat'),'MTparameters');
f_lMT_vMT_dM.save(fullfile(OutPath,'f_lMT_vMT_dM'));
f_ActivationDynamics.save(fullfile(OutPath,'f_ActivationDynamics'));
f_BackActivationDynamics.save(fullfile(OutPath,'f_BackActivationDynamics'));
f_FiberVelocity_TendonForce.save(fullfile(OutPath,'f_FiberVelocity_TendonForce'));
f_FiberLength_TendonForce.save(fullfile(OutPath,'f_FiberLength_TendonForce'));
f_forceEquilibrium_FtildeState.save(fullfile(OutPath,'f_forceEquilibrium_FtildeState'));
f_sumsqr_exp.save(fullfile(OutPath,'f_sumsqr_exp'));
f_T5.save(fullfile(OutPath,'f_T5'));
f_T4.save(fullfile(OutPath,'f_T4'));
f_T3.save(fullfile(OutPath,'f_T3'));




end

