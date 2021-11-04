function [] = Visualise2DModel(modelPath,motFile,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nRep = 1;
if ~isempty(varargin)
    nRep = varargin{1};
end

import org.opensim.modeling.*
m = Model(modelPath);
m.setUseVisualizer(true);
s = m.initSystem();
vis = m.getVisualizer();
vis.addDirToGeometrySearchPaths(char('C:\GBW_MyPrograms\OpenSim43\Geometry'))

dat = ReadMotFile(motFile);
dt = dat.data(2,1)-dat.data(1,1);
N = length(dat.data(:,1));

CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r',...
    'hip_flexion_l','lumbar_extension','knee_angle_r','knee_angle_l','ankle_angle_r','ankle_angle_l','mtp_angle_r','mtp_angle_l'};

% column index in datastrcture for every coordName
for i=1:length(CoordNamesAPI)
    IndexCoord(i) = find(strcmp(CoordNamesAPI{i},dat.names))-1;
end

% offset on x Coordinate
dat.data(:,2) = dat.data(:,2)-dat.data(1,2)-2;

% Qvect = s.updQ();
vis.show(s);
%% visualise the model
simbodyVis = vis.updSimbodyVisualizer();
simbodyVis.setShowSimTime(true);
for rep=1:nRep
    Qvect = s.getQ();
    for i=1:N
        dSel = dat.data(i,2:end);
        for j=1:length(CoordNamesAPI)
            if j==2 || j==3
                Qvect.set(j-1,dSel(IndexCoord(j)));
            else
                Qvect.set(j-1,dSel(IndexCoord(j))*pi/180);
            end
        end
        s.setQ(Qvect);
        pause(dt);
        vis.show(s);
    end
end



end

