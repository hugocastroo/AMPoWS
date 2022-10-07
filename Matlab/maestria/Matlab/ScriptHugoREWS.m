% Script to run OpenFAST with lidar simulator TEST HUGO


clearvars; close all; clc;

%% Load results from exe
resultData      = importdata('1p2_maininput.RO.dbg2');
% signalNames     = strsplit(resultData.textdata{end-1,1});
% signalNames     = signalNames(~cellfun(@isempty,signalNames));
signalNames     = resultData.textdata(2,:);
Time            = resultData.data(:,matches(signalNames,'Time'));
genPwr         = resultData.data(:,matches(signalNames,'VS_GenPwr'));
M_g             = resultData.data(:,matches(signalNames,'GenTq'));
genSpeed           = resultData.data(:,matches(signalNames,'GenSpeed'));
Omega           = resultData.data(:,matches(signalNames,'RotSpeed'));
WE_Vw           = resultData.data(:,matches(signalNames,'WE_Vw'));               %! Estimated wind speed [m/s]
WE_Vw_F         = resultData.data(:,matches(signalNames,'WE_Vw_F'));              %! Filtered estimated wind speed [m/s]
WE_VwI          = resultData.data(:,matches(signalNames,'WE_VwI'));              %! Integrated wind speed quantity for estimation [m/s]
WE_VwIdot       = resultData.data(:,matches(signalNames,'WE_VwIdot'));           %! Differentiated integrated wind speed quantity for estimation [m/s]
HorWindV        =  resultData.data(:,matches(signalNames,'HorWindV'));          %Horizontal wind Velocity

resultData1      = importdata('1p2_maininput.out');
signalNames1     = strsplit(resultData1.textdata{end-1,1});
signalNames1     = signalNames1(~cellfun(@isempty,signalNames1));
RtAeroCp         = resultData1.data(:,matches(signalNames1,'RtAeroCp'));

%% Plot results
%figure
controlador = 'RoscoOriginal';

ax1 = subplot(2,1,1);
hold on; grid on; box on
plot(Time(:,1),  WE_Vw,  '-','DisplayName','Estimated wind speed');
plot(Time(:,1), WE_Vw_F,  '-','DisplayName','Filtered estimated wind speed');
xlabel("Time in s");
ylabel("Velocity in m/s");
title("EWS-FEWS No-Delay");
lgd = legend;

ax2 = subplot(2,1,2);
hold on; grid on; box on
plot(Time(:,1),  WE_Vw,  '-','DisplayName','Estimated wind speed');
plot(Time(:,1)-2, WE_Vw_F,  '-','DisplayName','Filtered estimated wind speed');
xlabel("Time in s");
ylabel("Velocity in m/s");
title("EWS-FEWS Delayed");
lgd = legend;
linkaxes([ax1 ax2],'x');
xlim([0 600]);

figure

ax3 = subplot(2,1,1);
hold on; grid on; box on
plot(Time(:,1),  WE_Vw,  '-','DisplayName','Estimated wind speed');
plot(Time(:,1) , HorWindV,  '-','DisplayName','Hub height wind speed');
xlabel("Time in s");
ylabel("Velocity in m/s");
title("EWS - HorWindV No-Offset");
lgd = legend;

ax4 = subplot(2,1,2);
hold on; grid on; box on
coeff24hMA = ones(1, 160)/160;
HorWindVFiltered = filter(coeff24hMA, 1, HorWindV);
fDelay = (length(coeff24hMA)-1)/2;
plot(Time(:,1), WE_Vw,Time(:,1)-fDelay/160,HorWindVFiltered);
xlabel("Time in s");
ylabel("Velocity in m/s");
title("EWS - HorWindV Delayed");
legend('Estimated wind speed','Filtered Hub height wind speed');

linkaxes([ax3 ax4],'x');
xlim([0 600]);

%%
ro = 1.225;
radio = 64.5;
REWS = real(((2*genPwr)./(ro*pi*(radio^2)*RtAeroCp)).^(1/3));
coeff24hMA = ones(1, 160)/160;
fDelay = (length(coeff24hMA)-1)/2;
FREWS = filter(coeff24hMA, 1, REWS);
                                                  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   

% ax1 = subplot(1,1,1);                                                                                                                                                                                                                                                                            
% hold on; grid on; box on
% plot(Time(:,1),  genPwr,  '-','DisplayName','GenPwr');
% plot(Time(:,1), M_g.*genSpeed - 250000,  '-','DisplayName','GenPwrCalculated');
% xlabel("Time in s");
% ylabel("Power in N-m");
% title("Power");
% lgd = legend;
figure
ax5 = subplot(2,1,1);
hold on; grid on; box on
plot(Time(:,1), WE_Vw,Time(:,1),REWS);
xlabel("Time in s");
ylabel("Velocity in m/s");
title("EWS - Calculated REWS ");
legend('Estimated wind speed','Rotor Effective Wind Speed');

ax5 = subplot(2,1,2);
hold on; grid on; box on
plot(Time(:,1), WE_Vw,Time(:,1)-2*fDelay/160,FREWS);
xlabel("Time in s");
ylabel("Velocity in m/s");
title("EWS - FREWS ");
legend('Estimated wind speed','Filtered Rotor Effective Wind Speed');

xlim([0 600]);