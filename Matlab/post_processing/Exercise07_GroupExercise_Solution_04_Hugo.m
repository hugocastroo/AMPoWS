% -----------------------------
% Script: Generates a wind field with TurbSim and analizes it.
% Exercise 07 of Master Course 
% "Controller Design for Wind Turbines and Wind Farms"
% ------------
% Task:
% 
% ------------
% History:
% v01:  David Schlipf on 22-Nov-2021
% ----------------------------------
clearvars;clc;%close all;

%% b) Read the wind field into Matlab with readBLgrid.m.
nSeed = 6;
Seed_vec = 1:nSeed;
for URef = 4:2:24 %Loop for the wind speeds
    if(URef < 10)
        resolution = '%01d';
    else
        resolution = '%02d';
    end
    for iSeed = 1:nSeed
        Seed                = iSeed;
        TurbSimResultFile  	= ['e:\Tesis\Simulationen\Teil1\wind\NTM_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_turbsim.wnd'];
        %TurbSimResultFile  	= ['c:\Hugo\Temp\WorkingFolder\MeanSpectra\TurbulentWind\',num2str(Seed,'%03d'),'.wnd'];
        [velocity{iSeed}, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(TurbSimResultFile);
    end


%% c) Compare the analytic spectrum of a single point with an estimated one for the wind at hub height using pwelch.

% frequency
T       = size(velocity{1},1)*dt;
f_max   = 1/2*1/dt;
f_min   = 1/T;
df      = f_min;
f       = [f_min:df:f_max];
n_f     = length(f);

% SummVars
URef    = SummVars(3);
sigma   = SummVars(4)/100*URef;
h       = SummVars(1);

% Length scale IEC
L       = 8.1*42;

% time
t       = 0:dt:T-dt;
n_t     = length(t); 

% spectrum
S       = sigma^2 * 4*L/URef ./ (1 + 6 * f * L/URef ).^(5/3);

% modified spectrum to have the right std
S_mod   = S/(sum(S)*df)*sigma^2;

% estimation
n_Blocks                 	= 1;
SamplingFrequency           = 1/dt;
n_FFT                       = n_t/n_Blocks;
MyWindow                    = ones(n_FFT,1);
% MyWindow                    = hamming(n_FFT);

% idx for grid point
idx_y 	= 17;
idx_z  	= 17;

% plot
PlotSingleSpectra           = false;
figure
hold on;grid on;box on
set(gca,'xScale','log')
set(gca,'yScale','log')
xlabel('frequency [Hz]','FontSize', 20)
ylabel('Spectrum [(m/s)^2/Hz]','FontSize', 20)

for iSeed = 1:nSeed
    % extract signal
    u                           = velocity{iSeed}(:,1,idx_y,idx_z);

    % estimate spectrum
    [S_est(iSeed,:),f_est]   	= pwelch(u-mean(u),MyWindow,[],n_FFT,SamplingFrequency);
    
    % plot single spectrum
    if PlotSingleSpectra
        plot(f_est,S_est(iSeed,:),'.-')
    end
end

% calculate mean spectrum
S_est_mean          = mean(S_est);

% plot analytic and mean
plot(f,S,'r','Linewidth',1.5)
plot(f,S_mod,'b','Linewidth',1.5)
plot(f_est,S_est_mean,'k','Linewidth',1.5)

% legend etc.
% LegendCell          = {};
% if PlotSingleSpectra
%     LegendCell = cellstr(strcat('estimateSeed',num2str(Seed_vec','%03d')));
% end
% LegendCell{end+1}   = 'analytic';
% LegendCell{end+1}   = 'analytic+modified';
% LegendCell{end+1}   = 'mean';
% legend(LegendCell)
%title(['Wind speed ', num2str(URef,resolution), 'm/s ',num2str(idx_y),'Grid',num2str(nSeed),'Seeds']);
ylim([1e-2 1e3])
xlim([1e-3 1e1])
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
pdfpath = ['Windspeed_', num2str(URef,resolution), 'ms','_Gridcenter_',num2str(idx_y),'_Seeds',num2str(nSeed)];
saveas(gcf,pdfpath,'pdf');
end
