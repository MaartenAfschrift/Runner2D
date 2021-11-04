% This script provides an inital guess for the design variables. The guess is
% data-informed (DI). We use experimental data to provide an initial guess of
% the joint variables but set constant values to the muscle and back variables.
% We use a pre-defined guess for the final time.
%
% Author: Antoine Falisse
% Date: 9/9/2019
%
% Adapted to provide a guess for a full gait cycle as well
% Maarten Afschrift

function guess = getGuess_2D_DI_opti_mtp(d,Qs,nq,N,time_IC,NMuscle,jointi,scaling,varargin)

boolStride = 0;
if ~isempty(varargin)
    boolStride = varargin{1};
end

%% determine if we want to intial guess of a step or a full stride
if boolStride ==1
    Nin = N;
    N = Nin./2;
end

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs.allfilt));
Qs_spline.data(:,1) = Qs.allfilt(:,1);
Qdots_spline.data = zeros(size(Qs.allfilt));
Qdots_spline.data(:,1) = Qs.allfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs.allfilt));
Qdotdots_spline.data(:,1) = Qs.allfilt(:,1);
for i = 2:size(Qs.allfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allfilt(:,1),Qs.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allfilt(:,1),1);
end

% time vector
time    = Qs.allfilt(:,1);
nfr     = length(time);

%% Qs: data-informed
% Pelvis tilt
guess.Qs_all.data(:,jointi.pelvis.tilt) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis_tx
guess.Qs_all.data(:,jointi.pelvis.tx) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qs_all.data(:,jointi.pelvis.ty) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Hip flexion
guess.Qs_all.data(:,jointi.hip.l) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qs_all.data(:,jointi.hip.r) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Knee
guess.Qs_all.data(:,jointi.knee.l) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qs_all.data(:,jointi.knee.r) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qs_all.data(:,jointi.ankle.l) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qs_all.data(:,jointi.ankle.r) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% mtp joint (zeros)
guess.Qs_all.data(:,jointi.mtp.l) = zeros(nfr,1);
guess.Qs_all.data(:,jointi.mtp.r) = zeros(nfr,1);
% Trunk extension
guess.Qs_all.data(:,jointi.trunk.ext) = ...
    Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Interpolation
Qs_time = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'time'));
time_expi.Qs(1) = find(round(Qs_time,3) == round(time_IC(1),3));
time_expi.Qs(2) = find(round(Qs_time,3) == round(time_IC(2),3));
step = (Qs_time(time_expi.Qs(2))-Qs_time(time_expi.Qs(1)))/(N-1);
interval = Qs_time(time_expi.Qs(1)):step:Qs_time(time_expi.Qs(2));
guess.Qs = interp1(round(Qs_time,4),guess.Qs_all.data,round(interval,4));
guess.Qs(:,jointi.pelvis.tx) = guess.Qs(:,jointi.pelvis.tx) - ....
    guess.Qs(1,jointi.pelvis.tx);

%% Qdots: data-informed
% Pelvis tilt
guess.Qdots_all.data(:,jointi.pelvis.tilt) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis_tx
guess.Qdots_all.data(:,jointi.pelvis.tx) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdots_all.data(:,jointi.pelvis.ty) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Hip flexion
guess.Qdots_all.data(:,jointi.hip.l) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdots_all.data(:,jointi.hip.r) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Knee
guess.Qdots_all.data(:,jointi.knee.l) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdots_all.data(:,jointi.knee.r) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qdots_all.data(:,jointi.ankle.l) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdots_all.data(:,jointi.ankle.r) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% mtp joint
guess.Qdots_all.data(:,jointi.mtp.l) = zeros(nfr,1);
guess.Qdots_all.data(:,jointi.mtp.r) = zeros(nfr,1);
% Trunk extension
guess.Qdots_all.data(:,jointi.trunk.ext) = ...
    Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Interpolation
guess.Qdots = interp1(round(Qs_time,4),guess.Qdots_all.data,round(interval,4));

% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots: data-informed
% Pelvis tilt
guess.Qdotdots_all.data(:,jointi.pelvis.tilt) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis_tx
guess.Qdotdots_all.data(:,jointi.pelvis.tx) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdotdots_all.data(:,jointi.pelvis.ty) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Hip flexion
guess.Qdotdots_all.data(:,jointi.hip.l) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdotdots_all.data(:,jointi.hip.r) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Knee
guess.Qdotdots_all.data(:,jointi.knee.l) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdotdots_all.data(:,jointi.knee.r) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle
guess.Qdotdots_all.data(:,jointi.ankle.l) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdotdots_all.data(:,jointi.ankle.r) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% mtp angle
guess.Qdotdots_all.data(:,jointi.mtp.l) = zeros(nfr,1);
guess.Qdotdots_all.data(:,jointi.mtp.r) = zeros(nfr,1);
% Trunk extension
guess.Qdotdots_all.data(:,jointi.trunk.ext) = ...
    Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Interpolation
guess.Qdotdots = interp1(round(Qs_time,4),guess.Qdotdots_all.data,...
    round(interval,4));

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Arm variables
guess.a_b = 0.1*ones(N,nq.trunk);
guess.e_b = 0.1*ones(N,nq.trunk);

%% Final time
guess.tf = Qs.time(end);

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));

%% Add symmetric second step when selectin full stride guess
if boolStride ==1
    
    % add second stride for naive guess     
    guess.a         = [guess.a; guess.a];
    guess.FTtilde   = [guess.FTtilde; guess.FTtilde];
    guess.vA        = [guess.vA; guess.vA];
    guess.dFTtilde  = [guess.dFTtilde; guess.dFTtilde];
    guess.a_b       = [guess.a_b ; guess.a_b ];
    guess.e_b       = [guess.e_b; guess.e_b];
    guess.tf        = guess.tf*2;
    
    % LR symmetry for kinematics
    %1) 'pelvis_tilt',          (1) 
    %3)'pelvis_tx',             (2)
    %5)'pelvis_ty',...          (3)
    %7)'hip_flexion_l',         (4)
    %9)'hip_flexion_r'          (5)
    %11)'knee_angle_l'          (6)
    %13 'knee_angle_r'          (7)
    %15)'ankle_angle_l',        (8)
    %17)'ankle_angle_r',        (9)
    %19)'lumbar_extension'};    (10)
    %is1 = [1:6 7:18 19:20];
    isqdd = [1:3 5 4 7 6 9 8 10];   % flip L and R angles
    isq = [1:6 9 10 7 8 13 14 11 12 17 18 15 16 21:22 19:20 23:24]; % flip L and R joint angles
    guess.QsQdots = [guess.QsQdots; guess.QsQdots(:,isq)];
    guess.Qdotdots = [guess.Qdotdots; guess.Qdotdots(:,isqdd)];
    
    
    % adjust pelvis tx
    guess.QsQdots(N+1:end,jointi.pelvis.tx*2-1) = guess.QsQdots(N+1:end,jointi.pelvis.tx*2-1) + guess.QsQdots(N,jointi.pelvis.tx*2-1);
    % only rotation around z-axis for floatin base, so no adjusments in
    % roations needed here for the floating base in this 2D example
    N = Nin;
end


%% Update guess information for use with opti implementation
% (also evaluate at collocation points

% add values for last mesh point (states)
guess.QsQdots   = [guess.QsQdots; guess.QsQdots(end,:)];
guess.a         = [guess.a; guess.a(end,:)];
guess.FTtilde   = [guess.FTtilde; guess.FTtilde(end,:)];
guess.a_b       = [guess.a_b; guess.a_b(end,1)];

% approximate values at the collocation points
gmesh = linspace(0,1,N+1);
gcoll = linspace(0,1,N*d);
gmeshSlack = linspace(0,1,N);
guess.QsQdots_c = interp1(gmesh , guess.QsQdots,gcoll);
guess.a_c       = interp1(gmesh, guess.a,gcoll);
guess.FTtilde_c = interp1(gmesh, guess.FTtilde,gcoll);
guess.dFTtilde_c = interp1(gmeshSlack, guess.dFTtilde,gcoll);
guess.Qdotdots_c = interp1(gmeshSlack, guess.Qdotdots,gcoll);
guess.a_b_c = interp1(gmesh, guess.a_b,gcoll);

end