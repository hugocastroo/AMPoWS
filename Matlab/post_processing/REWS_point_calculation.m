%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a single point in a wind field.
%In this case, the point in the wind field, is the point in the center of
%wind field (33x33), the point in front o the hub. The calculation
%ilustrates the wind speed at Hub height since the point is in front of the
%hub. The REWS and the spectra are calculated and compared with the
%analytic values.
%------------------------------------------------------------
%V1.0 2022.10.28
%------------------------------------------------------------
close all; clear all; clc;
%%
%Add needed directories for the simulations:
addpath('WindFiles') %directory for the wind fields
%Create the wnd field using turbsim if needed
%CreateWindFields(4,8); %4 ms with seeds from 1 to 8
cd WindFiles\
dos(['TurbSim_x64.exe NTM_RandSeed1-401_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-402_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-403_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-404_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-405_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-406_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-407_turbsim.inp' ]);
dos(['TurbSim_x64.exe NTM_RandSeed1-408_turbsim.inp' ]);

%Read the wnd field from the turbsim wnd file and store the variables in
%Mat files
Load_wnd_files(4,8);
%% Calculate the analytic and estimate spectrum of a single point at hub height.
seeds = 8;
for URef = 4
    for i = 1:1:seeds %Loop for the seeds
        %Load the matlab variables
        [velocity, y, z, nz, ny, dz, dy, dt, z1, SummVars] = Load_mat_variables(URef,i);
        %Set the coordinates for the grid point
        idx_y               = 17;
        idx_z               = 17;
        %Read the velocity values for the grid point
        u                   = velocity(:,1,idx_y,idx_z);
        %Allocate variables for faster processing
        uTotal = zeros(length(u),seeds);
        uX = zeros(length(u),1);
        uTotal(:,i) = u; %Store the different speeds for the different seeds
        % frequency data
        T       = length(velocity)*dt;
        f_max   = 1/2*1/dt;
        f_min   = 1/T;
        df      = f_min;
        f       = [f_min:df:f_max];
        n_f     = length(f);
        
        % SummVars Get information from the file variables
        sigma   = SummVars(4)/100*URef;
        
        % Length scale IEC
        L       = 8.1*42;
        
        % time
        t       = 0:dt:T-dt;
        n_t     = length(t); 
        
        % spectrum calculation
        S       = sigma^2 * 4*L/URef ./ (1 + 6 * f * L/URef ).^(5/3);
        
        % estimation "More blocks mean better resolution, but less frequency lenght since the time series is beeing splitted into more parts"
        nBlocks                     = 8;
        SamplingFrequency           = 1/dt;
        n_FFT                       = n_t/nBlocks;
        [S_est,f_est]               = pwelch(u-mean(u),hamming(n_FFT),[],n_FFT,SamplingFrequency);
        
        % plot
        %figure
        hold on;grid on;box on
        plot(f_est,S_est,'.-')
        plot(f,S,'-')
        set(gca,'xScale','log')
        set(gca,'yScale','log')
        xlabel('frequency [Hz]')
        ylabel('Spectrum [(m/s)^2/Hz]')
        legend(['One point Hub estimate', ['velocity', int2str(URef)]],'Hub analytic')
    end
    %% Try making it with the mean of the 8 different signals
    for i=1:1:length(uTotal) %Loop for the seeds
    uX(i) = mean(uTotal(i,:));
    end
% % Try, concatnating the 8 different signals in just one
%     uX = [uTotal(:)];
%     nBlocks                     = 8;
%     SamplingFrequency           = 1/dt;
%     T       = length(uX)*dt;
%     t       = 0:dt:T-dt;
%     n_t = length(t);
%     n_FFT                       = n_t/nBlocks;
    %%
    [S_estTotal,f_estTotal]               = pwelch(uX-mean(uX),hamming(n_FFT),[],n_FFT,SamplingFrequency);
    %%figure
    hold on;grid on;box on
    plot(f_estTotal,S_estTotal,'.-')
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    xlabel('frequency [Hz]')
    ylabel('Spectrum [(m/s)^2/Hz]')
    legend(['One point Hub estimate', ['velocity', int2str(URef)]],'Hub analytic','Total')
end
