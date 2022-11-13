%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the an estimation of the REWS using the p
%welch Matlab function, after calculating the n estimations, it means the
%estimations to have a value closer to the analytic spectra. The more seeds
%the better
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------    
function [S_est_mean,f_est,WindFieldTurbSim,EstimationParam] = MeanPwelchEst(velocity,AnSpecParam, TurbsimParam, nSeed)

    %Rotos disc parameters
    R                   = 63;
    [Y,Z]               = meshgrid(TurbsimParam.y,TurbsimParam.z-AnSpecParam.h);
    DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
    IsInRotorDisc       = DistanceToHub<=R;
    % estimation parameters
    EstimationParam.nBlocks                     = 1;
    EstimationParam.SamplingFrequency           = 1/TurbsimParam.dt;
    EstimationParam.n_FFT                       = AnSpecParam.n_t/EstimationParam.nBlocks;
    EstimationParam.MyWindow                    = ones(EstimationParam.n_FFT,1);
    % MyWindow                    = hamming(n_FFT);
    
    %Create an array with the different spectrum estimations with pwelch using every different wind feld according to the different seeds
    S_est = zeros(nSeed,((EstimationParam.n_FFT/2)+1));                                     % allocation                                            
    WindFieldTurbSim = zeros(nSeed,EstimationParam.n_FFT);
    for iSeed = 1:nSeed
        v_0             = NaN(AnSpecParam.n_t,1);
        for i_t = 1:1:AnSpecParam.n_t
            CurrentWind     = squeeze(velocity{iSeed}(i_t,1,:,:));          % extract signal from wind field array and squeze it
            WindField       = CurrentWind(IsInRotorDisc);
            v_0(i_t,1)      = mean(WindField);
        end
        WindFieldTurbSim(iSeed,:) = v_0;
        % estimate spectrum
        [S_est(iSeed,:),f_est]   	= pwelch(v_0-mean(v_0),EstimationParam.MyWindow,[],EstimationParam.n_FFT,EstimationParam.SamplingFrequency);
    end

    %Calculate the mean REWS using the array with the different spectrum data from every seed
    S_est_mean          = mean(S_est);

return;