%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the an estimation of the REWS using the p
%welch Matlab function, after calculating the n estimations, it means the
%estimations to have a value closer to the analytic spectra. The more seeds
%the better
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------    
function [Scross_mean, FCross] = MeanCrossEst(EstimationParam, REWS, TEREWS, nSeed)

    %Cross Spectra calculation for n wind fields
    SCross_est = zeros(nSeed,((EstimationParam.n_FFT/2)+1));
    for iSeed = 1:nSeed
        % estimate cross spectrum with cpsd for n wind fields
        [SCross_est(iSeed,:),FCross] = cpsd(REWS(iSeed,:)-mean(REWS(iSeed,:)),TEREWS(iSeed,:)-mean(TEREWS(iSeed,:)),EstimationParam.MyWindow,[],EstimationParam.n_FFT,EstimationParam.SamplingFrequency);
    end
    
    %Get the mean spectrum from n seeds
    Scross_mean = mean(SCross_est);

return