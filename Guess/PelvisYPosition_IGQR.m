
%% Debug pelvis Y position
import casadi.*;

% external function
cd('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo\Sim_2D\ExternalFunctions');
S.ExternalFunc2 = 'PredSim_2D_Contact3_pp.dll';    % this one is with the pinjoint mtp 
F  = external('F',S.ExternalFunc2);
cd('C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo\Sim_2D\Guess');
% indexes
GRFi = [13:16];
jointi.pelvis.ty    = 3;

PyVect = 0.86:0.001:0.98;
OutVect = nan(16,length(PyVect));
for i=1:length(PyVect)
    Input = [zeros(12,1); zeros(12,1) ;zeros(12,1)];
    Input(jointi.pelvis.ty*2-1) = PyVect(i);
    Out = F(Input);
    OutVect(:,i) = full(Out);
end

figure();
FyPelvis = OutVect(jointi.pelvis.ty,:);
plot(PyVect,FyPelvis);
xlabel('Pelvis position');
ylabel('Residual Fy');
PelvisPos = PyVect(FyPelvis>=0);
disp(['Optimal pelvis position for this subjects is ' num2str(PelvisPos(1))]);
