%WIND_CONFIG_NTM Creates specific labels for the Normal-Turbulance-Model simulation 
%   (according to IEC61400) from basic parameter columns in configuration 
%   sheet by extending dlc_cell with necessary columns.
%
% Copyright (c) 2021 Hannah Dentzien, Ove Hagge Ellhöft
% Copyright (c) 2021 Jens Geisler

function [dlc_cell, turbsim_trig] = wind_config_NTM(dlc_cell,row) 


wind_speed = dlc_cell{row,find_label_or_create(dlc_cell,'Wind-Speed',true)}; 
duration = dlc_cell{row,find_label_or_create(dlc_cell,'Duration',true)};
seed = dlc_cell{row,find_label_or_create(dlc_cell,'Seed',true)};
turbclass = dlc_cell{row,find_label_or_create(dlc_cell,'Turb-Class',true)};
HubHt = dlc_cell{row,find_label_or_create(dlc_cell,'{Hub_Height}',true)};
slope = dlc_cell{row,find_label_or_create(dlc_cell,'Wind-Slope',true)};

% search for URef label (windspeed in TurbSim-Inputfile)
[idx,dlc_cell] = find_label_or_create(dlc_cell,'URef',false) ;
dlc_cell{row,idx}=wind_speed;  % write windspeed

% search for WindType label (Type of inputfile for inflowwind.dat-file)
[idx,dlc_cell] = find_label_or_create(dlc_cell,'WindType',false) ;
dlc_cell{row,idx}='3'; 

[idx,dlc_cell] = find_label_or_create(dlc_cell,'IECturbc',false) ;
dlc_cell{row,idx} = turbclass; 

% search for IEC_WindType label (TurbModel Turbsim-Inputfile)
[idx,dlc_cell] = find_label_or_create(dlc_cell,'IEC_WindType',false) ;
dlc_cell{row,idx}='NTM' ; 

% search for UsableTime (sim-time TurbSim-Inputfile)
%[idx,dlc_cell] = find_label_or_create(dlc_cell,'UsableTime',false) ;
%dlc_cell{row,idx}=duration; 

% search for AnalysisTime (sim-time TurbSim-Inputfile)
[idx,dlc_cell] = find_label_or_create(dlc_cell,'AnalysisTime',false) ;
dlc_cell{row,idx}=duration; 

% search for TMax label (sim-time main input (.fst))
[idx,dlc_cell] = find_label_or_create(dlc_cell,'TMax',false) ;
dlc_cell{row,idx}=duration;

% search for RandSeed1 label (First Random Seed Turbsim-Inputfile)
[idx,dlc_cell] = find_label_or_create(dlc_cell,'RandSeed1',false) ;
dlc_cell{row,idx} = seed;

[idx,dlc_cell] = find_label_or_create(dlc_cell,'HubHt',false) ;
dlc_cell{row,idx} = HubHt;

[idx,dlc_cell] = find_label_or_create(dlc_cell,'VFlowAng',false) ;
dlc_cell{row,idx} = slope;


turbsim_trig = true;

