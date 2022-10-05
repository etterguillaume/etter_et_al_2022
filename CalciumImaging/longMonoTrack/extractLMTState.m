function [LMT_state] = extractLMTState(optosignal);
%MSEXTRACTLMTSTATE Summary of this function goes here
%   Detailed explanation goes here

%% Params
threshold = 0.2;
Fs = 30;


%% Init vectors 
binary_signal = 0*optosignal;

binary_signal(optosignal>threshold)=1;
onset_signal = diff(binary_signal);
onset_signal(end+1) = 0;
flash_ts = find(onset_signal==1);
flash_int = diff(flash_ts/30);
flash_int(end+1)=0;
flash_freq = 1./flash_int;

blank_to_slow=[];
slow_to_fast=[];
fast_to_solid=[];

% First stim is obvious from deriv data
blank_to_slow(1) = 1;

for i = 2:length(flash_freq)
   if flash_freq(i-1) < 4 & flash_freq(i) > 4
       blank_to_slow(end+1)=i;
   elseif flash_freq(i-1) < 6 & flash_freq(i) > 6
       slow_to_fast(end+1)=i;
   elseif flash_freq(i-1) > 6 & flash_freq(i) < 4
       fast_to_solid(end+1) = i;
   end
end

%% Assumed drop in frequency corresponding to solid light
flash_ts(end+1)=flash_ts(end)+1
fast_to_solid(end+1) = length(flash_ts);

state = binary_signal*3+1;
state2write = 1;

frame_buffer = 10;

ct=1;

for i=1:length(state)-frame_buffer    
    if state2write == 1 && ct<=length(blank_to_slow) && flash_ts(blank_to_slow(ct))<i
        blank_to_slow(1)=[];
    end
    if state2write == 2 && ct<=length(slow_to_fast) && flash_ts(slow_to_fast(ct))<i
        slow_to_fast(1)=[];
    end
    if state2write == 3 && ct<=length(fast_to_solid) && flash_ts(fast_to_solid(ct))<i
        fast_to_solid(1)=[];
    end
        
    %% Re-order state-transitions
    while length(slow_to_fast)>=ct+frame_buffer && slow_to_fast(ct)<=blank_to_slow(ct)+frame_buffer
        slow_to_fast(1)=[];
    end
    while length(fast_to_solid)>=ct+frame_buffer && fast_to_solid(ct)<=blank_to_slow(ct)+frame_buffer
        fast_to_solid(1)=[];
    end
    while length(fast_to_solid)>=ct+frame_buffer && fast_to_solid(ct)<=slow_to_fast(ct)+frame_buffer
        fast_to_solid(1)=[];
    end
   
    %% Detect state transitions
   if state2write==1 && ct<=length(blank_to_slow) && flash_ts(blank_to_slow(ct))==i
      state2write = 2;
   elseif state2write==2 && ct<=length(slow_to_fast) && flash_ts(slow_to_fast(ct))==i
      state2write = 3;
   elseif state2write==3 && ct<=length(fast_to_solid) && flash_ts(fast_to_solid(ct))==i
      state2write = 4;
   elseif state2write==4 && state(i)==1
      state2write = 1;
      if ct+1<= length(blank_to_slow) & ct+1<= length(slow_to_fast) & ct+1<= length(fast_to_solid)
        ct=ct+1;
        while flash_ts(blank_to_slow(ct)) < i
        blank_to_slow(1)=[];
        end
      end
   end
       
   state(i)=state2write; 
end

LMT_state = state;

end

