%% Create New OsimModel


import org.opensim.modeling.*


%% Model with shorter hamstrings

m = Model('Gait18_vOsim43.osim');
mSel = m.getMuscles.get('hamstrings_r');
mSel.setTendonSlackLength(0.3);
mSel = m.getMuscles.get('hamstrings_l');
mSel.setTendonSlackLength(0.3);
m.print('Gait18_shorterHamstring.osim');