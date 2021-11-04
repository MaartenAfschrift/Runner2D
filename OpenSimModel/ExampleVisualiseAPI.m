% Example visualise model
clear all;
import org.opensim.modeling.*
m = Model('Gait18_Antoine2.osim');
m.setUseVisualizer(true);
s = m.initSystem();
vis = m.getVisualizer();
vis.addDirToGeometrySearchPaths(char('C:\GBW_MyPrograms\OpenSim43\Geometry'))
vis.show(s);


datafile = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D\Results\Example_Batch\DefaultSim_2D_speed_11_q.mot';
dat = ReadMotFile(datafile);
nfr = dat.data(:,1);
dt = dat.data(2,1)-dat.data(1,1);
N = length(dat.data(:,1));

% Coords = m.getCoordinateSet;
% for i=1:Coords.getSize()
%     CoordNamesAPI{i} = char(Coords.get(i-1).getName());
% end
CoordNamesAPI = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_r',...
    'hip_flexion_l','lumbar_extension','knee_angle_r','knee_angle_l','ankle_angle_r','ankle_angle_l','mtp_angle_r','mtp_angle_l'};

% column index in datastrcture for every coordName
for i=1:length(CoordNamesAPI)
    IndexCoord(i) = find(strcmp(CoordNamesAPI{i},dat.names))-1;
end

% Qvect = s.updQ();
% vis.show(s);
%%
% Qvect = s.getQ();
% Qvect.set(8,0.2); % 3 is heup R 4 is heup L, 5 is lumbar 6 is knie 7 knie L
% s.setQ(Qvect);
% vis.show(s);
simbodyVis = vis.updSimbodyVisualizer();
simbodyVis.setShowSimTime(true);
for rep=1:1
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
        s.setTime(dat.data(i,1));
        s.setQ(Qvect);
        pause(dt);
        vis.show(s);
%         simbodyVis.zoomCameraToShowAllGeometry()
    end
end
simbodyVis.pointCameraAt
% simbodyVis.delete
% simbodyVis.setShowFrameNumber(true)

