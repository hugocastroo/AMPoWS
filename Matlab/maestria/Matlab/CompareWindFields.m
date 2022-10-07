%% 1. Initialization
clear all; close all;clc;
%run('..\AddingWitlisPaths.m')

% 18 m/s seed 1801
PathToWND = 'c:\ChineseOpenFastTutorial\AMPoWS\examples\generated\wind\NTM_RandSeed1-1801_turbsim';
[velocity, ~, ~, ~, ~, dz, dy, dt, ~, ~, SummVars] = readBLgrid([PathToWND,'.wnd']);
[windfield1]         = velocity2windfield(velocity,dz,dy,dt,SummVars);

 % 3 m/s seed 201
% PathToWND = 'c:\ChineseOpenFastTutorial\AMPoWS\examples\generated\wind\NTM_RandSeed1-1201_turbsim';
% [velocity, ~, ~, ~, ~, dz, dy, dt, ~, ~, SummVars] = readBLgrid([PathToWND,'.wnd']);
% [windfield2]         = velocity2windfield(velocity,dz,dy,dt,SummVars);

% compare both
yIdx    = 33;
zIdx    = 33;
t       = windfield1.grid.t;
u1      = windfield1.u(yIdx,:,zIdx);
%u2      = windfield2.u(yIdx,:,zIdx);
%%
figure
hold on;box on;grid on
plot(t,u1,'DisplayName','NTMURef-20turbsim')
%plot(t,u2,'DisplayName','NTMURef-18turbsim')

lgd = legend;

figure
hold on;box on;grid on
plot(t,(u1-mean(u1))./std(u1,1),'DisplayName','NTMURef-20turbsim')
%plot(t,(u2-mean(u2))./std(u2,1),'DisplayName','NTMURef-18turbsim')

lgd = legend;
