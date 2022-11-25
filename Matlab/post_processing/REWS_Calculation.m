%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a wind field using n different seeds for n different wind fields.
%The wind field, has a grid size of 33x33 and is generated with TurbSim,
%sampling frequency 20Hz. The Estimated Turbsim REWS-Spectra is then
%compared with the Estimated REWS-Spectra from the ROSCO controller using a
%Kalman filter, sampling frequency 200Hz.
%------------------------------------------------------------
%V1.2 2022.11.19 - Time Delay function was added, Plot function was modified
%V1.1 2022.11.17 - Cross Spectra correction
%V1.0 2022.11.09 - Based on Script from David Schlip's lecture
% ----------------------------------
clearvars;clc;%close all;

% Set the basic info for the simulation
nSeed = 99;          %n Different seeds for each wind field speed
URefMinimum = 14;
URefMaximum = 24;
URefStep = 2;
tic
%% Loop for every different wind speed
for URef = URefMinimum:URefStep:URefMaximum
    fprintf(['Geschwindigkeit: ',num2str(URef)]);
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
    TEREWS = MatchFrequency(DataRosco,EstimationParam, nSeed);
 
    %Calculate the mean REWS using an overlapped signal
    [S_mean_TEREWS,f_est_TEREWS] = MeanPwelchEstRosco(AnSpecParam,TurbsimParam,TEREWS,nSeed);

    %Cross Spectra calculation for n wind fields
    [Scross_mean, f_Cross] = MeanCrossEst(EstimationParam, REWS, TEREWS, nSeed);

    %Plot Spectra data
    FrequenciesToPlot = {AnSpecParam.f,f_mean,f_est_TEREWS};
    SignalsToPlot = {S_RR,S_mean_REWS,S_mean_TEREWS};
    SignalNames = {'Analytic','Estimated mean REWS','Estimated mean TEREWS','Mean Cross'};
    xText = 'frequency [Hz]'; %description for the x axis
    yText = 'Spectrum [(m/s)^2/Hz]'; %description for the y axis
    tText = ''; %description for the title
    logViewx = true; %If true the plot is x axis of the plot is set to log scale
    logViewy = true; %If true the plot is y axis of the plot is set to log scale
    SavePlot = true; %If true, the plot is saved as pdf in full horizontal page size
    PlotSpectraResults(FrequenciesToPlot,SignalsToPlot,SignalNames,SavePlot,URef,xText,yText,tText,logViewx,logViewy)

    % Magnitude squared coherence (gamma square)
    cohsqr = (abs(Scross_mean)).^2./(S_mean_REWS.*S_mean_TEREWS);
    
    %Plot coherence data
    PlotSpectraResults({f_Cross},{cohsqr},{'GammaSqrSpectra'},true,URef,'frequency [Hz]','coherence [-]','',true,false)
    
    TimeDelay(URef/URefStep-1,:) = DelayCalculation(EstimationParam,REWS,TEREWS,nSeed);

    %Compare the windfields if needed
    %TestWindField = sqrt((TEREWS.^2+REWS.^2)/2); %Signal more similar to the turb sim wind field for a more similar spectrum test
%     time = 0:1/EstimationParam.SamplingFrequency:AnSpecParam.T-1/EstimationParam.SamplingFrequency; %Time vector
%     for i = 1:nSeed
%         TimeToPlot = {time+TimeDelay{URef/URefStep-1}{i},time};        
%         %SignalsToPlot = {TEREWS(i,:),REWS(i,:),TestWindField(i,:)};
%         %SignalNames = {'TEREWS','REWS','Test Modified signal'};
%         SignalsToPlot = {TEREWS(i,:),REWS(i,:)};
%         SignalNames = {'TEREWS','REWS'};
%         PlotSpectraResults(TimeToPlot,SignalsToPlot,SignalNames,false,URef,'time[s]','wind speed [m/s]',[num2str(URef), 'm/s'],false,false)
%     end
clear DataRosco; clear TurbSimWind; clear TEREWS; clear REWS;
end
toc
%Plot the time delay series if needed,
% for iSeed = 1:length(TimeDelay)
%     time{iSeed} = (1:nSeed);
%     windspeed{iSeed} = (abs(TimeDelay(iSeed,1:nSeed)));
%     names{iSeed} = ([num2str(URefMinimum+(URefStep*(iSeed-1))),' m/s']);
%     meansDelay(iSeed) = (abs(mean(TimeDelay(iSeed,1:nSeed))));
% end
% PlotSpectraResults(time,windspeed,names,false,URefMinimum+(URefStep*(iSeed-1)),'Seed [-]', 'Wind speed [m/s]','Time Delays',false,false)
% yline(meansDelay,'--',names,'Linewidth',2,'HandleVisibility','off');