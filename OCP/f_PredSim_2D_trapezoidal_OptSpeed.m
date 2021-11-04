function [Results] = f_PredSim_2D_trapezoidal_OptSpeed(S)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



%% Adding the casadi path seems to be needed to run processes in batch
AddCasadiPaths();

%% Default settings
S = GetDefaultSettings(S);

%% User inputs (typical settings structure)
% settings for optimization
N           = S.N;          % number of mesh intervals
W           = S.W;          % weights optimization

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathRepo = S.pathRepo;
addpath(genpath(pathRepo));
% Loading external functions.
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
F  = external('F',fullfile(pathExternalFunctions,S.ExternalFunc));

%% Indices external function
% Indices of the elements in the external functions
% External function: F
% First, joint torques.
jointi = getJointi();

% some vectors for later use
jointi.all          = jointi.pelvis.tilt:jointi.trunk.ext; % all
jointi.gr_pelvis    = jointi.pelvis.tilt:jointi.pelvis.ty; % ground-pelvis
nq.all              = length(jointi.all); % all
nq.abs              = length(jointi.gr_pelvis); % ground-pelvis
nq.leg              = 3; % #joints needed for polynomials
nq.trunk            = 1; % trunk
% External function: F1 (post-processing purpose only)
% Ground reaction forces (GRFs)
GRFi.r              = 13:14;
GRFi.l              = 15:16;
GRFi.all            = [GRFi.r,GRFi.l];
% model mass
body_mass = S.mass;

%% collocation scheme
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[~,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters

% Muscles from one leg
muscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};
NMuscle = length(muscleNames)*2;

% Indices of the muscles actuating the different joints for later use
tl = load(fullfile(pathRepo,'Polynomials',S.PolyFolder,'muscle_spanning_joint_INFO.mat'));
[~,mai] = MomentArmIndices_2D(muscleNames,tl.muscle_spanning_joint_INFO);

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
OutPath = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_lMT_vMT_dM = Function.load(fullfile(OutPath,'f_lMT_vMT_dM'));
% fgetMetabolicEnergySmooth2004all = Function.load(fullfile(OutPath,'fgetMetabolicEnergySmooth2004all'));
f_BackActivationDynamics = Function.load(fullfile(OutPath,'f_BackActivationDynamics'));
f_FiberVelocity_TendonForce = Function.load(fullfile(OutPath,'f_FiberVelocity_TendonForce'));
f_FiberLength_TendonForce = Function.load(fullfile(OutPath,'f_FiberLength_TendonForce'));
f_forceEquilibrium_FtildeState = Function.load(fullfile(OutPath,'f_forceEquilibrium_FtildeState'));
f_sumsqr_exp = Function.load(fullfile(OutPath,'f_sumsqr_exp'));
f_T5 = Function.load(fullfile(OutPath,'f_T5'));
f_T4 = Function.load(fullfile(OutPath,'f_T4'));
f_T3 = Function.load(fullfile(OutPath,'f_T3'));
f_ActivationDynamics = Function.load(fullfile(OutPath,'f_ActivationDynamics'));

fgetMetabolicEnergySmooth2004all = Function.load(fullfile(pathCasADiFunctions,'EnergyModels','fgetMetabolicEnergySmooth2004all'));
fgetMetabolicEnergy_MargariaSmooth  = Function.load(fullfile(pathCasADiFunctions,'EnergyModels', 'fgetMetabolicEnergy_MargariaSmooth'));
load(fullfile(OutPath,'MTparameters.mat'),'MTparameters');

%% additional muscle parameters
pctst = getSlowTwitchRatios_2D(muscleNames);
tension     = getSpecificTensions_2D(muscleNames);
pctst= [pctst; pctst];
[massM] =GetMuscleMass(MTparameters,tension');
massM = [massM massM];
Fmax = [MTparameters(1,:) MTparameters(1,:)];

%% Experimental data

% We extract experimental data to set bounds and initial guesses if needed
pathData = [pathRepo,'/OpenSimModel/',S.IKFolder];
joints = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_l',...
    'hip_flexion_r','knee_angle_l','knee_angle_r','ankle_angle_l',...
    'ankle_angle_r','lumbar_extension'};

%% Bounds
[bounds,scaling] = getBounds_2D_mtp_OptSpeed(NMuscle,nq,jointi);

%% Initial guess

% Extract joint positions from average walking motion
motion              = 'running';
nametrial.id        = ['motion_ave_',motion];
nametrial.IK        = ['IK_',nametrial.id];
pathIK              = [pathData,'/IK/',nametrial.IK,'.mot'];
Qs_IK               = getIK(pathIK,joints);
time_IC             = [Qs_IK.time(1),Qs_IK.time(end)];

% The initial guess depends on the settings
if S.IDig<3
    if S.IDig == 2
        guess = getGuess_2D_DI_opti_mtp(3,Qs_IK,nq,N,time_IC,NMuscle,jointi,scaling,S.stride);
    else
        guess = getGuess_2D_QR_opti_mtp(N,3,nq,NMuscle,scaling,S.v_tgt,jointi,S.IG_PelvisY);
    end
else
    % get guess based on previous simulation results
    GuessFolder     = fullfile(pathRepo,'Results',S.ResultsF_ig);
    pathIK          = fullfile(GuessFolder,[S.savename_ig '.mot']);
    Qs_ig           = getIK(pathIK,joints);
    if S.stride == 1
        Qs_ig.allfilt   = Qs_ig.allfilt(1:end,:);
        Qs_ig.time      = Qs_ig.time(1:end,:);
        time_IC         = [Qs_ig.time(1),Qs_ig.time(end)];
    else
        nfr = length(Qs_ig.allfilt);
        iend = nfr/2;
        Qs_ig.allfilt   = Qs_ig.allfilt(1:iend,:);
        Qs_ig.time      = Qs_ig.time(1:iend,:);
        time_IC         = [Qs_ig.time(1),Qs_ig.time(iend)];
    end
    % get the initial guess
    guess = getGuess_2D_DI_opti_mtp(3,Qs_ig,nq,N,time_IC,NMuscle,jointi,scaling,S.stride);
end

%% formulate the optimization problem

% create the opti structure
opti    = casadi.Opti();

% variables at mesh points
tf      = opti.variable(1,1);           % final time
a       = opti.variable(NMuscle,N+1);   % Muscle activations
FTtilde = opti.variable(NMuscle,N+1);   % Muscle-tendon forces
x       = opti.variable(2*nq.all,N+1);  % Qs and Qdots
ab      = opti.variable(1,N+1);         % back activations

% controls at mesh points
dFTtilde = opti.variable(NMuscle,N+1);    % time derivative muscle force
qdd      = opti.variable(nq.all,N+1);     % joint accelerations
e        = opti.variable(NMuscle,N+1);    % muscle excitations
eb       = opti.variable(1,N+1);          % back excitations

% bounds on variables at mesh points
opti.subject_to(bounds.tf.lower < tf < bounds.tf.upper);
opti.subject_to(bounds.a.lower'*ones(1,N+1)  < a  < bounds.a.upper'*ones(1,N+1));
opti.subject_to(bounds.FTtilde.lower'*ones(1,N+1)  < FTtilde  < bounds.FTtilde.upper'*ones(1,N+1));
opti.subject_to(bounds.QsQdots.lower'*ones(1,N+1) < x < bounds.QsQdots.upper'*ones(1,N+1));
opti.subject_to(bounds.a_b.lower  < ab  < bounds.a_b.upper);

% bound on initial state
opti.subject_to(x(3,1) == 0);   % pelvis Tx should be zero

% bounds on controls at mesh points
opti.subject_to(bounds.Qdotdots.lower'*ones(1,N+1) < qdd < bounds.Qdotdots.upper'*ones(1,N+1));
opti.subject_to(bounds.dFTtilde.lower'*ones(1,N+1) < dFTtilde < bounds.dFTtilde.upper'*ones(1,N+1));
opti.subject_to(bounds.e.lower < e < bounds.e.upper);
opti.subject_to(bounds.e_b.lower < eb < bounds.e_b.upper);

% Unscale the states and controlss
SFTtilde    = FTtilde.*(scaling.FTtilde'*ones(1,N+1));
Sx          = x.*(scaling.QsQdots'*ones(1,N+1));
SdFTtilde	= dFTtilde.*scaling.dFTtilde;
Sqdd        = qdd.*(scaling.Qdotdots'*ones(1,N+1));

% Provide expression for the distance traveled
pelvis_tx0 = Sx(2*jointi.pelvis.tx-1,1);
pelvis_txf = Sx(2*jointi.pelvis.tx-1,N+1);
dist_trav_tot = pelvis_txf-pelvis_tx0; % distance traveled

% pre-allocate matrices to store temporary variables (easy to get output)
% Out_lMT = MX(NMuscle,N);    Out_vMT = MX(NMuscle,N);
% Out_Tk  = MX(nq.all,N);     Out_GRF = MX(4,N);

% Time step
h = tf/N;

% quadrature function
J = 0;

% Loop over the mesh points
for k = 1:N
    % Get muscle-tendon lengths, velocities, and moment arms
    % Left leg
    qin_l = [Sx(jointi.hip.l*2-1,k),...
        Sx(jointi.knee.l*2-1,k),Sx(jointi.ankle.l*2-1,k)];
    qdotin_l = [Sx(jointi.hip.l*2,k),...
        Sx(jointi.knee.l*2,k),Sx(jointi.ankle.l*2,k)];
    [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM(qin_l,qdotin_l);
    MA_hip_l    =  MA_l(mai(1).mus.l',1);
    MA_knee_l   =  MA_l(mai(2).mus.l',2);
    MA_ankle_l  =  MA_l(mai(3).mus.l',3);
    % Right leg
    qin_r = [Sx(jointi.hip.r*2-1,k),...
        Sx(jointi.knee.r*2-1,k),Sx(jointi.ankle.r*2-1,k)];
    qdotin_r = [Sx(jointi.hip.r*2,k),...
        Sx(jointi.knee.r*2,k), Sx(jointi.ankle.r*2,k)];
    [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM(qin_r,qdotin_r);
    % Here we take the indices from left since the vector is 1:NMuscle/2
    MA_hip_r    = MA_r(mai(1).mus.l',1);
    MA_knee_r   = MA_r(mai(2).mus.l',2);
    MA_ankle_r  = MA_r(mai(3).mus.l',3);
    % Both legs
    lMTk_lr     = [lMTk_l;lMTk_r];          
%     Out_lMT(:,k) = lMTk_lr;
    vMTk_lr     = [vMTk_l;vMTk_r];          
%     Out_vMT(:,k) = vMTk_lr;
    % Call external function (Skeleton dynamics using .dll file)
    [Res]       = F([Sx(:,k);Sqdd(:,k)]);
    Tk          = Res(jointi.all);
%     Out_Tk(:,k) = Tk;
%     Out_GRF(:,k)= Res(GRFi.all);
    % add feedforward and feedback componets
    % backward euler integration scheme
    qk      = Sx(1:2:end,k);      qdk = Sx(2:2:end,k);    qddk = Sqdd(:,k);
    qk1     = Sx(1:2:end,k+1);    qdk1 = Sx(2:2:end,k+1); qddk1 = Sqdd(:,k+1); 
    dadt    = f_ActivationDynamics(e(:,k),a(:,k));
    dadt1    = f_ActivationDynamics(e(:,k+1),a(:,k+1));
    dadt_b   = f_BackActivationDynamics(eb(k),ab(k));
    dadt_b1   = f_BackActivationDynamics(eb(k+1),ab(k+1));
    Xk  = [qk;      qdk;      SFTtilde(:,k);      a(:,k);      ab(k)];         % vector with states
    Xdk  = [qdk;     qddk;     SdFTtilde(:,k);     dadt;        dadt_b];  % vector with state derivatives
    Xk1  = [qk1;     qdk1;    SFTtilde(:,k+1);    a(:,k+1);    ab(k+1)]; % vector with states (dt+1)
    Xdk1 = [qdk1;     qddk1;   SdFTtilde(:,k+1);   dadt1;       dadt_b1];  % vector with state derivatives
    
    opti.subject_to(TrapezoidalIntegrator(Xk,Xk1,Xdk,Xdk1,h) == 0);
                    
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffk,FTk,Fcek,Fpassk,Fisok] = ...
        f_forceEquilibrium_FtildeState(a(:,k),SFTtilde(:,k),...
        SdFTtilde(:,k),lMTk_lr,vMTk_lr);    
    % Compute joint moments generated by muscles
    Ft_hip_l    = FTk(mai(1).mus.l',1);
    T_hip_l     = f_T4(MA_hip_l,Ft_hip_l);
    Ft_hip_r    = FTk(mai(1).mus.r',1);
    T_hip_r     = f_T4(MA_hip_r,Ft_hip_r);
    Ft_knee_l   = FTk(mai(2).mus.l',1);
    T_knee_l    = f_T5(MA_knee_l,Ft_knee_l);
    Ft_knee_r   = FTk(mai(2).mus.r',1);
    T_knee_r    = f_T5(MA_knee_r,Ft_knee_r);
    Ft_ankle_l  = FTk(mai(3).mus.l',1);
    T_ankle_l   = f_T3(MA_ankle_l,Ft_ankle_l);
    Ft_ankle_r  = FTk(mai(3).mus.r',1);
    T_ankle_r   = f_T3(MA_ankle_r,Ft_ankle_r);
    
    % compute joint moment of mtp joint
    q_mtp_l = Sx(jointi.mtp.l*2-1,k); qdot_mtp_l = Sx(jointi.mtp.l*2,k);
    q_mtp_r = Sx(jointi.mtp.r*2-1,k); qdot_mtp_r = Sx(jointi.mtp.r*2,k);
    T_mtp_l = -S.k_mtp*q_mtp_l -S.d_mtp*qdot_mtp_l;
    T_mtp_r = -S.k_mtp*q_mtp_r -S.d_mtp*qdot_mtp_r;
    
    % Constraints: Muscle-driven joint torques
    opti.subject_to(Tk(jointi.hip.l,1) - T_hip_l ==0);
    opti.subject_to(Tk(jointi.hip.r,1) - T_hip_r ==0);
    opti.subject_to(Tk(jointi.knee.l,1) - T_knee_l ==0);
    opti.subject_to(Tk(jointi.knee.r,1) - T_knee_r ==0);
    opti.subject_to(Tk(jointi.ankle.l,1) - (T_ankle_l) ==0);     % exo torque + muscle torque should equal ID torque
    opti.subject_to(Tk(jointi.ankle.r,1) - (T_ankle_r) ==0);     % exo torque + muscle torque should equal ID torque
    opti.subject_to(Tk(jointi.mtp.l,1) - T_mtp_l == 0);
    opti.subject_to(Tk(jointi.mtp.r,1) - T_mtp_r == 0);
    opti.subject_to(Tk(jointi.trunk.ext,1)./scaling.BackTau - ab(:,k) == 0);
    % Constraints: Null pelvis residuals
    opti.subject_to(Tk(jointi.gr_pelvis,1) == 0);
    % constraints contraction dynamics
    opti.subject_to(Hilldiffk == 0);
    
%     if W.E>0
%         % get the metabolic energy
%         % Get muscle fiber lengths
%         [~,lMtildek] = f_FiberLength_TendonForce(SFTtilde(:,k),lMTk_lr);
%         % Get muscle fiber velocities
%         [vMk,~] = f_FiberVelocity_TendonForce(SFTtilde(:,k),SdFTtilde(:,k),lMTk_lr,vMTk_lr);
%         % Get metabolic energy rate        
%         [e_tot,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(a(:,k),a(:,k),...
%             lMtildek,vMk,Fcek,Fpassk,massM,pctst,Fisok,Fmax,body_mass,10);   
%         
% %         e_tot = fgetMetabolicEnergy_MargariaSmooth(Fcek,vMk,0.1);
%     else
%         e_tot= zeros(NMuscle,1);
%     end
    
    % objective function
    J = J + 1/(dist_trav_tot)*(...
        W.A*(f_sumsqr_exp(a(:,k),W.exp_A))*h + ...
        W.back*(sumsqr(eb(k)))*h +...
        W.acc*(sumsqr(qdd(:,k)))*h + ...
        W.A*(f_sumsqr_exp(e(:,k),W.exp_A))*h + ...
        W.u*(sumsqr(dFTtilde(:,k)))*h);
    
end

% constraint for perdiodic motion
if S.stride
    opti.subject_to(FTtilde(:,1) == FTtilde(:,end));
    opti.subject_to(a(:,1) == a(:,end));
    opti.subject_to(ab(:,1) == ab(:,end));
    opti.subject_to(x([1:2 4:end],1) == x([1:2 4:end],end));    % not for x-direction pelvis (only position)
else
    % constraint for symmetric motion
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    QsInvA = [jointi.pelvis.tilt:2*jointi.pelvis.tilt,...
        2*jointi.pelvis.tx,2*jointi.pelvis.ty-1:2*jointi.trunk.ext]';
    QsInvB = [jointi.pelvis.tilt:2*jointi.pelvis.tilt,...
        2*jointi.pelvis.tx,2*jointi.pelvis.ty-1:2*jointi.pelvis.ty,...
        2*jointi.hip.r-1:2*jointi.hip.r,...
        2*jointi.hip.l-1:2*jointi.hip.l,...
        2*jointi.knee.r-1:2*jointi.knee.r,...
        2*jointi.knee.l-1:2*jointi.knee.l,...
        2*jointi.ankle.r-1:2*jointi.ankle.r,...
        2*jointi.ankle.l-1:2*jointi.ankle.l,...
        2*jointi.mtp.r-1:2*jointi.mtp.r,...
        2*jointi.mtp.l-1:2*jointi.mtp.l,...
        2*jointi.trunk.ext-1,2*jointi.trunk.ext]';
    opti.subject_to(x(QsInvA,end) - x(QsInvB,1) == 0);
    opti.subject_to(a(:,end) - a(orderMusInv,1) == 0);
    opti.subject_to(ab(:,end) - ab(:,1) == 0);
    opti.subject_to(FTtilde(:,end) - FTtilde(orderMusInv,1) == 0);
end
% constraint on average speed
vel_aver_tot = dist_trav_tot/tf;
%J = 0.01.*J/N - 10*vel_aver_tot.^2;
J = 0.01.*J/N - 10*vel_aver_tot;
% opti.subject_to(vel_aver_tot - S.v_tgt == 0);

% set the initial guess
opti.set_initial(tf, guess.tf);
opti.set_initial(x, guess.QsQdots');
opti.set_initial(a, guess.a');
opti.set_initial(FTtilde, guess.FTtilde');
opti.set_initial(ab, guess.a_b');
opti.set_initial(e, 0.1);
opti.set_initial(dFTtilde,[guess.dFTtilde;guess.dFTtilde(end,:)]');
opti.set_initial(eb,[guess.e_b;guess.e_b(end,:)]');
opti.set_initial(qdd, [guess.Qdotdots;guess.Qdotdots(end,:)]');

% options for ipopt
optionssol.ipopt.hessian_approximation  = 'limited-memory';
optionssol.ipopt.mu_strategy            = 'adaptive';
optionssol.ipopt.linear_solver          = S.linear_solver;
optionssol.ipopt.max_iter               = 1500;
optionssol.ipopt.tol                    = 1*10^(-S.tol_ipopt);

% solve the OCP
opti.minimize(J);   % minimimze objective function
opti.solver('ipopt',optionssol);

% open diary
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
if ~isfolder(OutFolder)
    mkdir(OutFolder)
end
Outname = fullfile(OutFolder,[S.savename '_log.txt']);
diary(Outname);
% solve the OCP
% sol = opti.solve();
[w_opt,Results.stats] = solve_NLPSOL(opti,optionssol);
diary off

%% Unpack the optimal solution

ct = 1;
Results.tf = w_opt(ct);
ct = ct+1;

Results.a = reshape(w_opt(ct:ct+NMuscle*(N+1)-1),NMuscle,N+1)';
ct = ct + NMuscle*(N+1);

Results.F = reshape(w_opt(ct:ct+NMuscle*(N+1)-1),NMuscle,N+1)';
ct = ct + NMuscle*(N+1);

Results.x = reshape(w_opt(ct:ct+2*nq.all*(N+1)-1),2*nq.all,N+1)';
ct = ct + 2*nq.all*(N+1);

Results.ab = reshape(w_opt(ct:ct+1*(N+1)-1),1,N+1)';
ct = ct + 1*(N+1);

Results.dF = reshape(w_opt(ct:ct+NMuscle*(N+1)-1),NMuscle,N+1)';
ct = ct + NMuscle*(N+1);

Results.qdd = reshape(w_opt(ct:ct+nq.all*(N+1)-1),nq.all,N+1)';
ct = ct + nq.all*(N+1);

Results.e = reshape(w_opt(ct:ct+NMuscle*(N+1)-1),NMuscle,N+1)';
ct = ct + NMuscle*(N+1);

Results.eb = reshape(w_opt(ct:ct+1*(N+1)-1),1,N+1)';
ct = ct + 1*(N+1)-1;

if ct ~=length(w_opt)
   error('problem in extracting optimal solution'); 
end

%% unscale the results
Results.F   = Results.F.*(ones(N+1,1)*scaling.FTtilde);
Results.x   = Results.x.*(ones(N+1,1)*scaling.QsQdots);
Results.dF	= Results.dF.*(ones(N+1,1)*scaling.dFTtilde);
Results.qdd = Results.qdd.*(ones(N+1,1)*scaling.Qdotdots);

%% get speed

% walking speed
pelvis_tx0 = Results.x(1,2*jointi.pelvis.tx-1);
pelvis_txf = Results.x(end,2*jointi.pelvis.tx-1);
dist_trav_tot = pelvis_txf-pelvis_tx0; % distance traveled
Results.speed = dist_trav_tot./Results.tf;

% Muscle states
% Get muscle-tendon lengths, velocities, and moment arms
qin_l = Results.x(:,[jointi.hip.l jointi.knee.l  jointi.ankle.l]*2-1);
qdotin_l = Results.x(:,[jointi.hip.l jointi.knee.l  jointi.ankle.l]*2);
qin_r = Results.x(:,[jointi.hip.r jointi.knee.r  jointi.ankle.r]*2-1);
qdotin_r = Results.x(:,[jointi.hip.r jointi.knee.r jointi.ankle.r]*2);
Results.lMT = zeros(N,NMuscle);   
Results.vMT = zeros(N,NMuscle);   
Results.MA = zeros(N,NMuscle,3);
Results.lMtilde = zeros(N,NMuscle);
Results.lM = zeros(N,NMuscle);
Results.vMtilde = zeros(N,NMuscle);
Results.vM = zeros(N,NMuscle);
for i=1:N    
    % muscle tendon lengths
    [lMT_lt,vMT_lt,MA_lt] = f_lMT_vMT_dM(qin_l(i,:),qdotin_l(i,:));
    [lMT_rt,vMT_rt,MA_rt] = f_lMT_vMT_dM(qin_r(i,:),qdotin_r(i,:));
    Results.lMT(i,:) = [full(lMT_lt); full(lMT_rt)]' ;
    Results.vMT(i,:) = [full(vMT_lt); full(vMT_rt)]';
    Results.MA(i,1:NMuscle/2,:) = full(MA_lt);
    Results.MA(i,NMuscle/2+1:NMuscle,:) = full(MA_rt);    
    % muscle states
    [lMk,lMtildek] = f_FiberLength_TendonForce(Results.F(i,:),Results.lMT(i,:));
    [vMk,vMtildek] = f_FiberVelocity_TendonForce(Results.F(i,:),Results.dF(i,:),...
        Results.lMT(i,:),Results.vMT(i,:));
    Results.lMtilde(i,:) = full(lMtildek);
    Results.lM(i,:) = full(lMk);
    Results.vMtilde(i,:) = full(vMtildek);
    Results.vM(i,:) = full(vMk);
end

% get the m
Results.muscleNames = muscleNames;


%% save the results file

save(fullfile(OutFolder,[S.savename '.mat']),'Results');

%% export the simulated motion (simple for now)

% For visualization in OpenSim GUI
% Mesh points
tgrid = linspace(0,Results.tf,N+1);
q    = Results.x(:,1:2:end);
qOut = q;
JointAngle.labels = {'time','pelvis_tilt','pelvis_tx','pelvis_ty',...
    'hip_flexion_l','hip_flexion_r','knee_angle_l','knee_angle_r',...
    'ankle_angle_l','ankle_angle_r','mtp_angle_l','mtp_angle_r','lumbar_extension'};
% convert from radians to degress (not for the x and y translation)
qOut(:,[1 4:12]) = qOut(:,[1 4:12])*180./pi;

% add multiple gait cycles
QsInvA_q = 4:11;
QsInvB_q = [5 4 7 6 9 8 11 10];
nhs = 2;
if nhs >1
    qcopy = qOut(1:N,:);
    qOut = repmat(qcopy,nhs,1);
    dt = tgrid(end);
    dx = q(end,2);
    tgrid = repmat(tgrid(1:N),1,nhs);
    for i=2:nhs
        iv = (N)*(i-1)+1:(N)*i;
        qOut(iv,2) = qOut(iv,2)+dx*(i-1);    % pelvis tilt
        tgrid(iv) = tgrid(iv)+ dt*(i-1);        % time vector
        if rem(i,2) == 0
            % change sides
            qOut(iv,QsInvA_q) = qcopy(:,QsInvB_q);
        end
    end
end
JointAngle.data = [tgrid' qOut];
MotFile = fullfile(OutFolder,[S.savename '_q.mot']);
write_motionFile(JointAngle, MotFile);

if S.BoolVideo
    nRep = 3;
    Visualise2DModel(S.ModelPath,MotFile,nRep);
end




end

