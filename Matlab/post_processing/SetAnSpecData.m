%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script set the frequency data and parameters needed to
%calculate one analytic spectrum, based on the velocity length of the wind
%field and the basic parameters from turbsim
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------

function AnSpecParam = SetAnSpecData(TurbSimWind, TurbSimParam, URef)

    % After reading the files create the analytic spectrum of the wind speed make an estimation for the whole wind field using pwelch.
    % frequency data
    AnSpecParam.T       = size(TurbSimWind{1},1)*TurbSimParam.dt;
    AnSpecParam.f_max   = 1/2*1/TurbSimParam.dt;
    AnSpecParam.f_min   = 1/AnSpecParam.T;
    AnSpecParam.df      = AnSpecParam.f_min;
    AnSpecParam.f       = AnSpecParam.f_min:AnSpecParam.df:AnSpecParam.f_max;
    AnSpecParam.n_f     = length(AnSpecParam.f);
    
    % SummVars from wnd files
    AnSpecParam.sigma   = TurbSimParam.SummVars(4)/100*URef;
    AnSpecParam.h       = TurbSimParam.SummVars(1);
    
    % Length scale IEC
    AnSpecParam.L       = 8.1*42;
    
    % time
    AnSpecParam.t       = 0:TurbSimParam.dt:AnSpecParam.T-TurbSimParam.dt;
    AnSpecParam.n_t     = length(AnSpecParam.t); 

return;