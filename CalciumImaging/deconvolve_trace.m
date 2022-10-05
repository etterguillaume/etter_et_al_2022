function deconvolved_trace = deconvolve_trace(ca_trace,ca_time)
%MSEXTRACTFIRING Summary of this function goes here
%   Detailed explanation goes here

decayTime = 1;

V.Ncells = 1;
V.T = length(ca_trace);
V.dt = mode(diff(ca_time));
V.fast_plot = 0;

    F = ca_trace';
    temp = min(F);
    F = F-temp;
    %P.sig = mean(mad(F,1)*1.4826); 
    %P.gam = (1-V.dt/decayTime)*ones(V.Ncells,1); %tau=dt/(1-gam)
    %P.lam =  .2*ones(V.Ncells,1); %expected spikers per frame
    P.a = median(F,2); %size of dF per spike
    %P.b = -1*temp*ones(V.Ncells,1);%quantile(F,0.05); %baseline dF

    [n_best p_best V2 C] = fast_oopsi(F,V,P);
    deconvolved_trace = n_best;
%     plot(F)
%     hold on
%     plot(n_best,'r')
%     hold off
%     ylim([-.1 .5])
%     drawnow
%     segNum

% ms.n = n_best;

end