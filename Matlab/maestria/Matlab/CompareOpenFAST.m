% Script to run OpenFAST with lidar simulator
% using the exe 
% based on https://github.com/andrew-platt/openfast/tree/f/LidarSim 
% David Schlipf and Florian Thomas on 22-Mar-2021
% (c) sowento GmbH 
% Update DS on 24-Aug-2022 for Hugo: 
% - use FB controller to see Torque init issue

%clearvars; close all; clc;

%% Run exe
dos('openfast_x64.exe NRELOffshrBsline5MW_EOG16.fst');

%% Load results from exe
resultData      = importdata('NRELOffshrBsline5MW_EOG16.out');
signalNames     = strsplit(resultData.textdata{end-1,1});
signalNames     = signalNames(~cellfun(@isempty,signalNames));
Time            = resultData.data(:,contains(signalNames,'Time'));
u_HH            = resultData.data(:,contains(signalNames,'Wind1VelX'));
x_T_dot         = resultData.data(:,contains(signalNames,'TTDspFA'));
M_g             = resultData.data(:,contains(signalNames,'GenTq'));
theta           = resultData.data(:,contains(signalNames,'BldPitch1'));
Omega           = resultData.data(:,contains(signalNames,'RotSpeed'));

%% Plot results
%figure
controlador = 'RoscoOriginal';

ax1 = subplot(5,1,1);
hold on; grid on; box on
plot(Time,  u_HH,  '-','DisplayName',controlador);
xlabel("Time in s");
ylabel("Velocity in m/s");
title("Wind speed at hub height");
lgd = legend;
ax2 = subplot(5,1,2);
hold on; grid on; box on
plot(Time,  theta,  '-','DisplayName',controlador);
xlabel("Time in s");
ylabel("angle in deg");
title("Collective pitch angle");
lgd = legend;

ax3 = subplot(5,1,3);
hold on; grid on; box on
plot(Time,  M_g,  '-','DisplayName',controlador);
grid on;
xlabel("Time in s");
ylabel("Torque in kNm");
title("Generator torque");
lgd = legend;

ax4 = subplot(5,1,4);
hold on; grid on; box on
plot(Time,  Omega,  '-','DisplayName',controlador);
grid on;
xlabel("Time in s");
ylabel("Speed in rpm");
title("Rotor speed");
lgd = legend;

ax5 = subplot(5,1,5);
hold on; grid on; box on
plot(Time,  x_T_dot,  '-','DisplayName',controlador);
grid on;
xlabel("Time in s");
ylabel("Velocity in m/s");
title("Tower top fore-aft speed");

lgd = legend;

linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
xlim([0 60]);
