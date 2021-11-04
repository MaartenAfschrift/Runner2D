function [massM] =GetMuscleMass(params,tension)
%UNTITLED5 Summary of this function goes here
% 
FMo = params(1,:);
lMo = params(2,:);
volM = FMo.*lMo;
massM = volM.*(1059.7)./(tension*1e6);

end

