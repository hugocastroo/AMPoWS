%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script matches the frequencies between a ROSCO REWS signal
%and a wind field signal from TurbSim
%------------------------------------------------------------
%V1.0 2022.11.17 - HC
% ----------------------------------    
function [TEREWS] = MatchFrequency(DataRosco,EstimationParam, nSeed)

    signalaux = zeros(nSeed,120002);
    signalaux1 = zeros(nSeed,120000);
    for iSeed = 1:nSeed
        %signalaux(iSeed,:) = [DataRosco{iSeed}.WE_Vw(120001:end-1); DataRosco{iSeed}.WE_Vw(20001:120000)];
        signalaux(iSeed,:) = circshift(DataRosco{iSeed}.WE_Vw(20000:end),20000);
        signalaux1(iSeed,:) = signalaux(iSeed,1:end-2);
%         Y = circshift(DataRosco{iSeed}.WE_Vw(20000:end-1),20000);
%         figure
%         hold on; grid on;
%         plot(1:length(Y),Y);
%         plot(1:length(DataRosco{iSeed}.WE_Vw),DataRosco{iSeed}.WE_Vw)
    end
    
    TEREWS = zeros(nSeed,EstimationParam.n_FFT);       % allocation
    for iSeed = 1:nSeed
        for i = 1:(length(TEREWS))
            TEREWS(iSeed,i) = signalaux1(iSeed,i*10-9);
        end
    end
return