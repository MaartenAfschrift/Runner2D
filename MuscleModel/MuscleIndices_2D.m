% This function returns the muscle indices used in the optimization problem
% as compared to all muscles of the gait1018 model
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function musi = MuscleIndices_2D(muscleNames)
   
muscleNames_all = {'hamstrings_r','bifemsh_r','glut_max_r',...
    'iliopsoas_r','rect_fem_r','vasti_r','gastroc_r','soleus_r',...
    'tib_ant_r'};
    
count = 1;
musi = zeros(1,length(muscleNames));
for i = 1:length(muscleNames)       
    if (find(strcmp(muscleNames_all,muscleNames{i})) ~= 0)        
        musi(count) = find(strcmp(muscleNames_all,muscleNames{i}));
        count = count + 1;
    end
end

end
