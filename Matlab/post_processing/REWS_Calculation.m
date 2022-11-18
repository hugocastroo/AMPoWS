%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a wind field using n different seeds for n different wind fields.
%The wind field, has a grid size of 33x33 and is generated with TurbSim,
%sampling frequency 20Hz. The Estimated Turbsim REWS-Spectra is then
%compared with the Estimated REWS-Spectra from the ROSCO controller using a
%Kalman filter, sampling frequency 200Hz.
%------------------------------------------------------------
%V1.1 2022.11.17 - Cross Spectra correction
%V1.0 2022.11.09 - Based on Script from David Schlip's lecture
% ----------------------------------
clearvars;clc;%close all;
%% Set the basic info for the simulation
nSeed = 6;          %n Different seeds for each wind field speed
URefMinimum = 4;
URefMaximum = 24;
URefStep = 2;
%%
%Loop for every different wind speed
for URef = URefMinimum:URefStep:URefMaximum

    % Read the wnd files from Turbsim
    [TurbSimWind, TurbsimParam] = ReadWndFiles(URef,nSeed);
    
    % Analytic Spectrum - frequency data calculation using the wind field length and turbsim parameters
    AnSpecParam = SetAnSpecData(TurbSimWind, TurbsimParam, URef);
        
    % Spectra rotor-effective wind speed - rotor model
    S_RR = SpectraREWSCalculation(AnSpecParam, TurbsimParam, URef);
    
    % Mean pwelch calculation over n seeds for the TurbSim Wind fields
    [S_mean_REWS,f_mean,REWS,EstimationParam] = MeanPwelchEst(TurbSimWind, AnSpecParam, TurbsimParam, nSeed);

    %Read the n ROSCO dbg file for URef wind speed and store it in a cell array
    DataRosco = ReadRoscoDbgFiles(URef,nSeed);

    % Array modification, since frequencies in Turbsim(20 Hz) and in OpenFAST(200 Hz) are different.
    % Overlapping the signal, since the last 100 seconds shouldbe overlapped, because the rosco signal is 700 seconds long
    TEREWS = MatchFrequency(DataRosco,EstimationParam, nSeed);

    %Compare the windfields if needed
%         %Signal more similar to the turb sim wind field for a more similar spectrum test
%     TestWindField = sqrt((TEREWS.^2+REWS.^2)/2);
%     figure
%     hold on;grid on;box on
%     plot(1:length(TEREWS),mean(TEREWS),'DisplayName','Rosco overlapped wind field')
%     plot(1:length(REWS),mean(REWS),'DisplayName','TurbSim wind field')
%     plot(1:length(REWS),(mean(TestWindField)),'DisplayName','Test Modified signal');
%     xlabel('time[s]','FontSize', 20)
%     ylabel('wind speed [m/s]','FontSize', 20)
%     lgd = legend;
 
    %Calculate the mean REWS using an overlapped signal
    [S_mean_TEREWS,f_est_TEREWS] = MeanPwelchEstRosco(AnSpecParam,TurbsimParam,TEREWS,nSeed);

    %Cross Spectra calculation for n wind fields
    [Scross_mean, f_Cross] = MeanCrossEst(EstimationParam, REWS, TEREWS, nSeed);

    %Plot data
    SignalNames = {'Analytic','Estimated mean REWS','Estimated mean TEREWS','Mean Cross'};
    SavePlot = false; %If true, the plot is saved as pdf in full horizontal page size
    FrequenciesToPlot = {AnSpecParam.f,f_mean,f_est_TEREWS,f_Cross};
    SignalsToPlot = {S_RR,S_mean_REWS,S_mean_TEREWS,abs(Scross_mean)};
    PlotSpectraResults(FrequenciesToPlot,SignalsToPlot,SignalNames,SavePlot,URef);

    %% Magnitude squared coherence (gamma square)
    cohsqr = (abs(Scross_mean)).^2./(S_mean_REWS.*S_mean_TEREWS);

    %cohsqrCorrection = cohsqr - ((1/100)*(1-cohsqr).^2);
    
    k=(2*pi*f_Cross)/URef;

    grl = abs((Scross_mean)./(S_mean_REWS));

    % Coherence
    Distance    = 1;
    kappa       = 12*((f_Cross/URef).^2+(0.12/AnSpecParam.L).^2).^0.5;
    gamma_uu    = exp(-kappa.*Distance); % coherence between point 1 and 2 in u

    %Plot Coherence
    figure
    hold on;grid on;box on
    plot(k,cohsqr,'DisplayName','GammaSqrSpectra');
    plot(k,gamma_uu.^2,'DisplayName','Analytic Kaimal');
    ylim([0 1])
    xlim([10^-3.1 10^1])
    set(gca,'xScale','log')
    xlabel('k [rad/m]','FontSize', 20)
    ylabel('coherence [-]','FontSize', 20)
    title([num2str(URef), 'm/s']);
    lgd = legend;

    W_cutoff = 0.033333411319;
    W_delay = 2*pi*0.1;
    W_n = W_delay/W_cutoff;
    T_filter = (arctan(W_n))/(W_delay);

    %Plot transfer function
%     figure
%     hold on;grid on;box on
%     plot(k,grl,'DisplayName','Grl');
%     ylim([0 1])
%     xlim([10^-3.1 10^1])
%     set(gca,'xScale','log')
%     xlabel('k [rad/m]','FontSize', 20)
%     ylabel('GRL [-]','FontSize', 20)
%     title([num2str(URef), 'm/s']);
%     lgd = legend;

end
