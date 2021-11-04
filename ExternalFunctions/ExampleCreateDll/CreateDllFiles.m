
% add path with matlab code to create .dll files
addpath('C:\Users\u0088756\Documents\FWO\Software\GitProjects\CreateDll_PredSim');

% path to cpp file
CppDir = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D\ExternalFunctions\ExampleCreateDll';
Names = dir(fullfile(CppDir,'*.cpp'));

% additional path information
OsimSource = 'C:\opensim-ad-core-source';
OsimBuild = 'C:\GBW_MyPrograms\opensim-ad-core-build2';
DllPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\PredSim_2D\ExternalFunctions';
ExtFuncs = 'C:\opensim-ExternalFunc';
VSinstall = 'C:\Program Files (x86)\Microsoft Visual Studio 14.0';

% number of input arguments for external function
nInputDll = 36;

% create the dll file
for i = 1:length(Names)
    Name = Names(i).name;
    Name = Name(1:end-4); % remove .cpp
    CreateDllFileFromCpp(CppDir,Name,OsimSource,OsimBuild,DllPath,ExtFuncs,VSinstall,nInputDll);
end

