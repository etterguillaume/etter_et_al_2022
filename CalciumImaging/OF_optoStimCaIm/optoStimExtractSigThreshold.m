function [OpticSigThreshold] = optoStimExtractSigThreshold(behav)
%OPTOSTIMEXTRACTSIGTHRESHOLD Summary of this function goes here
%   Detailed explanation goes here

optoSignal = behav.optosignal;
smoothedOptoSignal = smooth(optoSignal, 'sgolay');
normalizedOptoSignal = smoothedOptoSignal-min(smoothedOptoSignal);
normalizedOptoSignal = normalizedOptoSignal./max(normalizedOptoSignal);

figure; plot(normalizedOptoSignal);
[~,OpticSigThreshold] = ginput(1);

end

