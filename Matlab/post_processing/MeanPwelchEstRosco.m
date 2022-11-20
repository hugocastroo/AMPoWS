%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the an estimation of the REWS - Estimate
%from the KF in the ROSCO Controller
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------    
function [S_mean_TEREWS,f_est_TEREWS] = MeanPwelchEstRosco(AnSpecParam,TurbsimParam,TEREWS,nSeed)

    % estimation parameters
    nBlocks                     = 1;
    SamplingFrequency           = 1/TurbsimParam.dt;
    n_FFT                       = AnSpecParam.n_t/nBlocks;
    MyWindow                    = ones(n_FFT,1);
    % MyWindow                    = hamming(n_FFT);

    %Create an array with the different spectrum estimations with pwelch using every different wind feld according to the different seeds
    S_estRosco = zeros(nSeed,((n_FFT/2)+1));                                     % allocation                                            
    for iSeed = 1:nSeed
        vRosco = TEREWS(iSeed,:);
        % estimate spectrum
        [S_estRosco(iSeed,:),f_est_TEREWS]   	= pwelch(vRosco-mean(vRosco),MyWindow,[],n_FFT,SamplingFrequency);
    end

    %Calculate the mean REWS using the array with the different spectrum data from every seed
    S_mean_TEREWS          = mean(S_estRosco);

    %Plot just one signal
%      for iSeed = 1:nSeed
%          plot(f_estRosco,S_estRosco(iSeed,:),'Linewidth',0.1,'DisplayName',['ROSCO seed',num2str(iSeed)])                                %Plot the estimation made with the mean values and pwelch
%      end

return;