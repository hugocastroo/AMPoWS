%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script matches the frequencies between a ROSCO REWS signal
%and a wind field signal from TurbSim
%------------------------------------------------------------
%V1.0 2022.11.17 - HC
% ----------------------------------    
function [TEREWS] = MatchFrequency(DataRosco,EstimationParam, nSeed)

    Roscof = 1/DataRosco{1}.Time(2)-DataRosco{1}.Time(1);   %Rosco sampling frequency
    FASTf = round(EstimationParam.SamplingFrequency);       %FAST sampling frequency
    Begin = (floor(length(DataRosco{1}.WE_Vw)/(floor(length(DataRosco{1}.WE_Vw)/Roscof)/100)))+2; %Get the array begin for the new signal, according to the length of the original signal and the sampling frequency
    TEREWS = zeros(nSeed,EstimationParam.n_FFT);       % allocation

    for i = 1:nSeed
        % Array modification, since frequencies in Turbsim(20 Hz) and in OpenFAST(200 Hz) are different.
        TEREWS(i,:) = DataRosco{i}.WE_Vw(Begin:Roscof/FASTf:end);
    end

return

