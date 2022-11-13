%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a wind field using n different seeds for different wind fields.
%The wind field, has a grid size of 33x33. The REWS and the spectra are calculated and compared with the
%analytic values.
%------------------------------------------------------------
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

    %Read the wnd files from Turbsim
    [velocity, TurbsimParam] = ReadWndFiles(URef,nSeed);
    
    % Analytic Spectrum - frequency data calculation using the wind field length and turbsim parameters
    AnSpecParam = SetAnSpecData(velocity, TurbsimParam, URef);
        
    % Spectra rotor-effective wind speed - rotor model
    S_RR = SpectraREWSCalculation(AnSpecParam, TurbsimParam, URef);
    
    % Mean pwelch calculation over n seeds
    [S_est_mean,f_est,WindFieldTurbSim] = MeanPwelchEst(velocity,AnSpecParam, TurbsimParam, nSeed);

    %Read the n ROSCO dbg file for URef wind speed and store it in a cell array
    DataRosco = ReadRoscoDbgFiles(URef,nSeed);

    % Array modification, since frequencies in Turbsim(20 Hz) and in OpenFAST(200) are different.
    %Overlapping the signal, since the last 100 seconds should not be used,
    %because the rosco signal is 700 seconds long
    for iSeed = 1:nSeed
        signalaux(iSeed,:) = [DataRosco{iSeed}.WE_Vm(120001:end-1); DataRosco{iSeed}.WE_Vm(20001:120000)];
    end
      
    ModifiedRoscoOverlapped = zeros(nSeed,12000);       % allocation
    for iSeed = 1:nSeed
        for i = 1:(length(ModifiedRoscoOverlapped))
            ModifiedRoscoOverlapped(iSeed,i) = signalaux(iSeed,i*10-9);
        end
    end

    %Calculate the mean REWS using an overlapped signal - BETTER RESULTS
    [S_est_meanRosco,f_estRosco] = MeanPwelchEstRosco(AnSpecParam,TurbsimParam,ModifiedRoscoOverlapped,nSeed);
    
    %%Cross Spectra Parameters
    vRoscoTotal = mean(ModifiedRoscoOverlapped);
    v_0Total = mean(WindFieldTurbSim);
    nBlocks                     = 1;
    SamplingFrequency           = 1/TurbsimParam.dt;
    n_FFT                       = AnSpecParam.n_t/nBlocks;
    MyWindow                    = ones(n_FFT,1);
    %Cross Spectra calculation
    [SCross(1,:),FCross] = cpsd(v_0Total-mean(v_0Total),vRoscoTotal-mean(vRoscoTotal),MyWindow,[],n_FFT,SamplingFrequency);
    
    %Plot data
    SignalNames = {'Analytic','Estimated mean Turbsim','Estimated mean Rosco','Cross Spectra'};
    SavePlot = false; %If true, the plot is saved as pdf in full horizontal page size
    FrequenciesToPlot = {AnSpecParam.f,f_est,f_estRosco,FCross};
    SignalsToPlot = {S_RR,S_est_mean,S_est_meanRosco,abs(SCross)};
    PlotSpectraResults(FrequenciesToPlot,SignalsToPlot,SignalNames,SavePlot,URef);

    %%
    %%Coherence calculation to compare the mscohere funcion with the 
    [gamma_Sq_est,fcoh_est]= mscohere(v_0Total-mean(v_0Total),vRoscoTotal-mean(vRoscoTotal),MyWindow,[],n_FFT,SamplingFrequency);

    % Magnitude squared coherence (gamma square)
    cohsqr = (abs(SCross)).^2./(S_est_mean.*S_est_meanRosco);
    k=(2*pi*FCross)/URef;

    % Coherence
    Distance    = 1 ;
    kappa       = 12*((FCross/URef).^2+(0.12/AnSpecParam.L).^2).^0.5;
    gamma_uu    = exp(-kappa.*Distance); % coherence between point 1 and 2 in u

    figure
    hold on;grid on;box on
    plot(k,cohsqr,'DisplayName','GammaSqrSpectra');
    plot(fcoh_est,gamma_Sq_est,'DisplayName','GammaSqrEstmscohere');
    plot(k,gamma_uu.^2,'DisplayName','Analytic');
    ylim([0 1])
    set(gca,'xScale','log')
    xlabel('wave number [rad/m]')
    ylabel('coherence [-]')
    title([num2str(URef), 'm/s']);
    lgd = legend;

end

%% 

% /////////////////////////////////// calculating the estimate without
% overlapping the signals
% Array modification, since frequencies in Turbsim(20 Hz) and in OpenFAST(200) are different.
%     ModifiedRoscoWE_Vw = zeros(nSeed,floor((length(DataRosco{nSeed}.WE_Vw)-1)/11.666));       % allocation
%     for iSeed = 1:nSeed
%         for i = 1:(length(ModifiedRoscoWE_Vw))
%             ModifiedRoscoWE_Vw(iSeed,i) = DataRosco{iSeed}.WE_Vw(floor(i*10)-9);
%         end
%     end
% 
%         %Calculate the mean REWS using the array with the different spectrum data from every seed
%     [S_est_meanRosco,f_estRosco] = MeanPwelchEstRosco(AnSpecParam,TurbsimParam,ModifiedRoscoWE_Vw,nSeed);
