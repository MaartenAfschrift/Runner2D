function [] = FitPolynomials(MainPath,ModelName,Modelpath,PolyFolder,Bool_RunMA)
%FitPolynomials Main function to fit polynomials on muscle analysis data.
%   Detailed explanation goes here



%% Settings

% Import the opensim API
import org.opensim.modeling.*

% flow control
% BoolEvaluatePolynomials_Debug = false;

% Output folder for this subject
OutFolder = PolyFolder;
SubjFolder = fullfile(MainPath,'Polynomials',OutFolder);
MA_path=fullfile(SubjFolder,'MuscleAnalysis');
if ~isfolder(MA_path)
    mkdir(MA_path)
end

%% Create dummy motion and run muscle analysis
if Bool_RunMA
    
    % bounds on dofs to create training dataset
    Bound_hipflex = [-50 50];
    Bound_knee = [-90 0];
    Bound_ankle = [-30 30];
%     Bound_mtp = [-30 10];
    
    % get the model coordinates
    m = Model(Modelpath);
    CoordSet = m.getCoordinateSet();
    nc = CoordSet.getSize();
    NamesCoordinates = cell(1,nc);
    for i = 1:nc
        NamesCoordinates{i} = char(CoordSet.get(i-1).getName());
    end
    % construct a file with generalized coordinates
    headers=[{'time'} NamesCoordinates];
    
    if strcmp(ModelName,'Gait18')
        
        % create the dummy motion
        n=5000; p=3;
        X = lhsdesign(n,p);
        X_scale=[diff(Bound_hipflex) diff(Bound_knee) diff(Bound_ankle)];
        X_min=[Bound_hipflex(1) Bound_knee(1) Bound_ankle(1)];
        Angles=X.*(ones(n,1)*X_scale)+(ones(n,1)*X_min);
        IndexAngles = [4:6]+1; % +1 because of time vector
        
        % path with dummy motion
        time=(1:n)./100;
        data=zeros(length(time),length(headers));
        data(:,1)=time;
        data(:,IndexAngles) = Angles;   % right leg
        data(:,IndexAngles+4) = Angles; % left leg
        pathDummyMotion = fullfile(SubjFolder,'dummy_motion_mtp.mot');
        generateMotFile(data,headers,pathDummyMotion);
        
        % select coordinate names to run muscle analysis
        CoordNames = headers(IndexAngles);
    end
    %Run a muscle analysis on the dummy motion
    disp('Muscle analysis running....');
    OpenSim_Muscle_Analysis(pathDummyMotion,Modelpath,MA_path,[time(1) time(end)],CoordNames);
end


%% import the muscle analysis data
% subject pre-fix
SubjPre = 'dummy_motion_mtp';
lMT = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_Length.sto']);
MA.hip.flex = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_hip_flexion_r.sto']);
MA.knee.flex = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_knee_angle_r.sto']);
MA.ankle.flex = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_ankle_angle_r.sto']);
% MA.mtp = ReadMotFile([MA_path,'/',SubjPre '_MuscleAnalysis_MomentArm_mtp_angle_r.sto']);
pathDummyMotion = fullfile(SubjFolder,'dummy_motion_mtp.mot');
dummy_motion = ReadMotFile(pathDummyMotion);
% Order of dofs: hip flex r knee flex r, ankle flex
% r, hip flex l, hip add l, hip rot l, knee flex l, ankle flex l, lumbar
% ext, lumbar bend, lumbar rot, subtalar r, subtalar l, mtp_r, mtp_l
order_Qs = [5:7];
q = dummy_motion.data(:,order_Qs).*(pi/180);

MuscleData.dof_names = dummy_motion.names(order_Qs);
muscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};
MuscleData.muscle_names = muscleNames;
lMTMat = nan(length(q),length(muscleNames));
dMMat = nan(length(q),length(muscleNames),3);
for m = 1:length(muscleNames)
    lMTMat(:,m)     = lMT.data(:,strcmp(lMT.names,muscleNames{m}));            % lMT
    dMMat(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.names,muscleNames{m}));    % hip_flex
    dMMat(:,m,2)    = MA.knee.flex.data(:,strcmp(lMT.names,muscleNames{m}));   % knee
    dMMat(:,m,3)    = MA.ankle.flex.data(:,strcmp(lMT.names,muscleNames{m}));  % ankle
%     dMMat(:,m,4)    = MA.mtp.data(:,strcmp(lMT.names,muscleNames{m}));         % mtp
end
MuscleData.lMT = lMTMat;
MuscleData.dM = dMMat;
MuscleData.q = q;
MuscleData.qdot = zeros(size(q));

%% Call PolynomialFit
max_order = 5;
[muscle_spanning_joint_INFO,MuscleInfo] = ...
    PolynomialFit_mtp(MuscleData,max_order);
save(fullfile(SubjFolder,'MuscleData.mat'),'MuscleData')
save(fullfile(SubjFolder,'muscle_spanning_joint_INFO.mat'),'muscle_spanning_joint_INFO')
save(fullfile(SubjFolder,'MuscleInfo.mat'),'MuscleInfo');


%% evaluate polynomials
% if BoolEvaluatePolynomials_Debug
%     import casadi.*
%     NMuscle_pol = 9;
%     musi_pol = 1:9;
%     nq.leg = 3;
%     muscle_spanning_info_m = muscle_spanning_joint_INFO(musi_pol,:);
%     MuscleInfo_m.muscle    = MuscleInfo.muscle(musi_pol);
%     q     = SX.sym('q',1,nq.leg);
%     qd    = SX.sym('qd',1,nq.leg);
%     lMT   = SX(NMuscle_pol,1);
%     vMT   = SX(NMuscle_pol,1);
%     dM    = SX(NMuscle_pol,nq.leg);
%     for i=1:NMuscle_pol
%         index_dof_crossing  = find(muscle_spanning_info_m(i,:)==1);
%         order               = MuscleInfo_m.muscle(i).order;
%         [mat,diff_mat_q]    = n_art_mat_3_cas_SX(q(1,index_dof_crossing),order);
%         lMT(i,1)            = mat*MuscleInfo_m.muscle(i).coeff;
%         vMT(i,1)            = 0;
%         dM(i,1:nq.leg)      = 0;
%         nr_dof_crossing     = length(index_dof_crossing);
%         for dof_nr = 1:nr_dof_crossing
%             dM(i,index_dof_crossing(dof_nr)) = ...
%                 (-(diff_mat_q(:,dof_nr)))'*MuscleInfo_m.muscle(i).coeff;
%             vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*...
%                 qd(1,index_dof_crossing(dof_nr)));
%         end
%     end
%     f_lMT_vMT_dM = casadi.Function('f_lMT_vMT_dM',{q,qd},{lMT,vMT,dM},...
%         {'q','qd'},{'lMT','vMT','dM'});
%     
%     lMTsel = zeros(5000,NMuscle_pol);
%     dMsel = zeros(5000,NMuscle_pol,3);
%     for i=1:5000
%         [lMTtemp,  ~, dMtemp] = f_lMT_vMT_dM(MuscleData.q(i,:),MuscleData.qdot(i,:));
%         lMTsel(i,:) = full(lMTtemp);
%         dMsel(i,:,:) = full(dMtemp);
%     end
%     
%     
%     figure(); 
%     subplot(1,2,1)
%     plot(MuscleData.q(:,2),MuscleData.lMT(:,6),'or','MarkerSize',1); hold on;
%     plot(MuscleData.q(:,2),lMTsel(:,6),'ok','MarkerSize',1);
%     subplot(1,2,2)
%     plot(MuscleData.q(:,2),squeeze(MuscleData.dM(:,6,2)),'or','MarkerSize',1); hold on;
%     plot(MuscleData.q(:,2),squeeze(dMsel(:,6,2)),'ok','MarkerSize',1)
%     
%     
% end
% 
% m = 6;









