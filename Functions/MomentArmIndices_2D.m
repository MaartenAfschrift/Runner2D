% This function returns the indices of the muscles actuating the different 
% joints (for use with the moment arms) and the indices of the muscles in 
% the vector containing all muscle names.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [Indmusi,mai] = ...
    MomentArmIndices_2D(muscleNames,muscle_spanning_joint_INFO)

NMuscle = length(muscleNames)*2;
iJoints = muscle_spanning_joint_INFO;
iActuated =  find(sum( muscle_spanning_joint_INFO) > 0);
for i = 1:length(muscleNames)
    % Muscles are ordered as left first and right second
    Indmusi.(muscleNames{i}(1:end-2)).l = ...
        find(strcmp(muscleNames,muscleNames{i}));
    Indmusi.(muscleNames{i}(1:end-2)).r = ...
        Indmusi.(muscleNames{i}(1:end-2)).l + NMuscle/2;
    % We select the muscles that are actuating the joints
    % In muscle_spanning_joint_INFO: the columns are organized as follows:
    % 1) hip flex, 2) knee flex, 3) ankle flex 4) mtp
    for j = iActuated
        if (muscle_spanning_joint_INFO(i,j) == 1)
           mai(j).mus.l(1,i) = Indmusi.(muscleNames{i}(1:end-2)).l;
           mai(j).mus.r(1,i) = Indmusi.(muscleNames{i}(1:end-2)).r;
        end        
    end    
end

for j = iActuated
    mai(j).mus.l(mai(j).mus.l == 0) = [];
    mai(j).mus.r(mai(j).mus.r == 0) = [];
end

end
