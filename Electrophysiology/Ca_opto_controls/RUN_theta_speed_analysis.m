%% Load the data:
load('/speed_GRIN_controls.mat') 

%% First, extract velocities from DLC extracted position data:
PV1100_baseline_velocity = extract_DLC_velocity(PV1100_baseline);
PV1102_baseline_velocity = extract_DLC_velocity(PV1102_baseline);
PV1103_baseline_velocity = extract_DLC_velocity(PV1103_baseline);
PV1045_baseline_velocity = extract_DLC_velocity(PV1045_baseline);
PV1046_baseline_velocity = extract_DLC_velocity(PV1046_baseline);
PV239_baseline_velocity = extract_DLC_velocity(PVJ20_239_baseline);
PV241_baseline_velocity = extract_DLC_velocity(PVJ20_241_baseline);

%% Construct speed vector
bin_size = 1;
speed_bins = 2:bin_size:14;

%% Extract theta speed correspondance
PV1045_theta_speed = extract_DLC_velocity_theta_correspondance(data_PV1045_OF_8Hz_638nm_5ms_1mW_5sONOFF.LFP.Data(:,2),PV1045_baseline_velocity,speed_bins);
PV1046_theta_speed = extract_DLC_velocity_theta_correspondance(data_PV1046_OF_8Hz_638nm_5ms_1mW_5sONOFF.LFP.Data(:,6),PV1046_baseline_velocity,speed_bins);
PV239_theta_speed = extract_DLC_velocity_theta_correspondance(data_239_OF_habituation.LFP.Data(:,4),PV239_baseline_velocity,speed_bins);
PV241_theta_speed = extract_DLC_velocity_theta_correspondance(data_241_OF_habituation.LFP.Data(:,1),PV241_baseline_velocity,speed_bins);
PV1100_theta_speed = extract_DLC_velocity_theta_correspondance(data_PV1100_baseline_noscope.LFP.Data(:,3),PV1100_baseline_velocity,speed_bins);
PV1102_theta_speed = extract_DLC_velocity_theta_correspondance(data_PV1102_baseline_noscope.LFP.Data(:,2),PV1102_baseline_velocity,speed_bins);
PV1103_theta_speed = extract_DLC_velocity_theta_correspondance(data_PV1103_baseline_nominiscope.LFP.Data(:,4),PV1103_baseline_velocity,speed_bins);

plot(PV1045_theta_speed,'DisplayName','PV1045_theta_speed');hold on;plot(PV1046_theta_speed,'DisplayName','PV1046_theta_speed');plot(PV1100_theta_speed,'DisplayName','PV1100_theta_speed');plot(PV1102_theta_speed,'DisplayName','PV1102_theta_speed');plot(PV1103_theta_speed,'DisplayName','PV1103_theta_speed');plot(PV239_theta_speed,'DisplayName','PV239_theta_speed');plot(PV241_theta_speed,'DisplayName','PV241_theta_speed');hold off;