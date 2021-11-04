function [IC1i,HS1] = DetectHeelContactPredSim(GRF,threshold,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% increase threshold untill you have at least on frame above the threshold
nFramesBelow= sum(GRF(:,2)<threshold);
while nFramesBelow == 0
    threshold = threshold + 1;
    nFramesBelow= sum(GRF(:,2)<threshold);
end
if threshold <100
    phase_tran_tgridi = find(GRF(:,2)<threshold,1,'last');
else
    % heelstrike is in between left and right leg simulation
    if threshold >100 && GRF(end,5)<20
        phase_tran_tgridi = find(GRF(:,2)<threshold,1,'last');
    else
        % heelstrike is on the left leg
        phase_tran_tgridi =[];
        threshold = 20;
    end
end
if ~isempty(phase_tran_tgridi)
    if phase_tran_tgridi == N
        temp_idx = find(GRF(:,2)>threshold,1,'first');
        if ~isempty(temp_idx)
            if temp_idx-1 ~= 0 && ...
                    find(GRF(temp_idx-1,2)<threshold)
                phase_tran_tgridi_t = temp_idx;
                IC1i = phase_tran_tgridi_t;
                HS1 = 'r';
            end
        else
            IC1i = phase_tran_tgridi + 1;
            HS1 = 'r';
        end
    else
        IC1i = phase_tran_tgridi + 1;
        HS1 = 'r';
    end
end



end

