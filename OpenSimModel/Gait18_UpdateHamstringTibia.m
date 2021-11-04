%% Create New OsimModel


import org.opensim.modeling.*


%% Model with shorter hamstrings

m = Model('gait10dof18musc_OsimInstall.osim');
mSel = m.getMuscles.get('hamstrings_r');
mSel.setTendonSlackLength(0.3);
mSel = m.getMuscles.get('hamstrings_l');
mSel.setTendonSlackLength(0.3);

mSel = m.getMuscles.get('tib_ant_r');
mSel.setMaxIsometricForce(3000);
mSel = m.getMuscles.get('tib_ant_l');
mSel.setMaxIsometricForce(3000);
m.print('Gait18_UpdateHamstringTibia.osim');