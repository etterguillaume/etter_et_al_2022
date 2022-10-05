function [cell_firing_characteristics] = msExtractTrans(ms)
%MSEXTRACTTRANS Summary of this function goes here
%   Detailed explanation goes here
%% Attention: all time data are expressed in samples!

%% Parameters
dt = median(diff(ms.time))/1000; % Conversion from ms to s
Fs = 1/dt;
dFF_threshold = 0.5; % 0.05
max_rise_time = 5; % in s
max_decay_time = 10; % in s

[bFilt,aFilt] = butter(2,  2/(30/2), 'low');

for trace_i = 1:size(ms.RawTraces,2);
    
rise_time = [];
decay_time = [];
prominence = [];
width = [];
max_rise_slope = [];
max_decay_slope = [];
rise_area = [];
decay_area = [];
frequency = [];
   

   raw_trace = ms.RawTraces(:,trace_i);
    filt_trace = zscore(filtfilt(bFilt,aFilt,raw_trace));
    current_trace = filt_trace;
  
   
   %% Peak detection
   [peaks,locs, width, prominence] = findpeaks(current_trace,'MinPeakHeight', dFF_threshold);
   
   if ~isempty(peaks) %& median(diff(locs))*dt > 1;
       for peak_i = 1:length(peaks);

       %% Rise time information
       if locs(peak_i) > round(max_rise_time*Fs);
       rise_period = current_trace(round(locs(peak_i)-max_rise_time*Fs):locs(peak_i));
       else
       rise_period = current_trace(1:locs(peak_i));
       end
       
       rise_deriv=diff(rise_period);
       rise_deriv(end+1)=0;
       
       mirror_rise_deriv = flipud(rise_deriv);
       
       zero_crossing=[];
       
       for rise_i = 2:length(rise_period);
           if mirror_rise_deriv(rise_i) < 0 & mirror_rise_deriv(rise_i-1) >= 0;
           zero_crossing(length(zero_crossing)+1) = rise_i;
           end
       end

       if ~isempty(zero_crossing) & zero_crossing(1)~=length(rise_deriv);
       rise_time(peak_i) = zero_crossing(1);
       max_rise_slope(peak_i) = max(rise_deriv(length(rise_deriv)-rise_time(peak_i):end));
       rise_area(peak_i) = trapz(rise_deriv(length(rise_deriv)-rise_time(peak_i):end)); %Numerical integration using trapezoid function
       else   
       %display('The rise time of transients is larger than the window size.');
       rise_time(peak_i) = length(rise_period);
       max_rise_slope(peak_i) = max(rise_deriv);
       rise_area(peak_i) = trapz(rise_period); %Numerical integration using trapezoid function
       end

       %% Decay time information
       if locs(peak_i) > length(ms.time) - max_decay_time*Fs;
       decay_period = current_trace(locs(peak_i):end);
       else
       decay_period = current_trace(locs(peak_i):locs(peak_i)+round(max_decay_time*Fs));
       end

       decay_deriv=diff(decay_period);
       decay_deriv(end+1)=0;

       zero_crossing=[];
       for rise_i = 2:length(decay_period);
           if decay_deriv(rise_i) > 0 & decay_deriv(rise_i-1) < 0;
           zero_crossing(length(zero_crossing)+1) = rise_i;
           end
       end
       
       if ~isempty(zero_crossing);
       decay_time(peak_i) = zero_crossing(1);
       max_decay_slope(peak_i) = max(decay_deriv(1:zero_crossing(1)));
       decay_area(peak_i) = trapz(decay_period(1:zero_crossing(1)));
       else   
       %display('The decay time of transients is larger than the window size.');
       decay_time(peak_i) = length(decay_period);
       max_decay_slope(peak_i) = max(decay_deriv);
       decay_area(peak_i) = trapz(decay_period);
       end
       
       end
       
   
   cell_firing_characteristics{trace_i}.start_time = locs-rise_time';
   cell_firing_characteristics{trace_i}.start_time(cell_firing_characteristics{trace_i}.start_time==0)=1; %% to avoid 0 indexing
   cell_firing_characteristics{trace_i}.end_time = locs+decay_time';
   cell_firing_characteristics{trace_i}.peak_idx = locs;
   cell_firing_characteristics{trace_i}.rise_time = rise_time';
   cell_firing_characteristics{trace_i}.decay_time = decay_time';
   cell_firing_characteristics{trace_i}.prominence = prominence;
   cell_firing_characteristics{trace_i}.width = width;
   cell_firing_characteristics{trace_i}.rise_area = rise_area';
   cell_firing_characteristics{trace_i}.decay_area = decay_area';
   cell_firing_characteristics{trace_i}.max_rise_slope = max_rise_slope';
   cell_firing_characteristics{trace_i}.max_decay_slope = max_decay_slope';
   
   cell_firing_characteristics{trace_i}.mean_rise_time = mean(rise_time');
   cell_firing_characteristics{trace_i}.mean_decay_time = mean(decay_time');
   cell_firing_characteristics{trace_i}.mean_prominence = mean(prominence);
   cell_firing_characteristics{trace_i}.mean_width = mean(width);
   cell_firing_characteristics{trace_i}.mean_rise_area = mean(rise_area');
   cell_firing_characteristics{trace_i}.mean_decay_area = mean(decay_area');
   
   else

   cell_firing_characteristics{trace_i} = [];
   
   end
   
   cell_firing_characteristics{trace_i}.frequency = mean(1./diff(cell_firing_characteristics{trace_i}.start_time));

   
end
%cell_firing_characteristics.detection_threshold = dFF_threshold;


end

