%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: Calculation of the time delay for every wind field using the xcorr function, the mean top
%values are then averaged and stored in the n+1 element of every array
%------------------------------------------------------------
%V1.0 2022.11.19 - HC
% ----------------------------------   
function TimeDelay = DelayCalculation(EstimationParam,REWS,TEREWS,nSeed)
    
    TimeDelay = zeros(1,nSeed); %allocate
    for i = 1:nSeed
        signal1= detrend(REWS(i,:),'constant'); %REWS Windfield
        signal2= detrend(TEREWS(i,:),'constant'); %TERES Windfield
        [c,lags]=xcorr(signal1,signal2,'coeff');
        top = max(c); %Top value in the Y axis
        TimeDelay(i) = interp1(c,lags*1/EstimationParam.SamplingFrequency,top); %Interpolation to get the X value (time) according to the Y axis
    end
    TimeDelay(nSeed+1) = mean((TimeDelay));

return