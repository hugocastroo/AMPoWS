%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the spectra of a rotor effective wind speed
%(rotor model) using the analytic spectra parameters and turbsim parameters
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------    
function S_RR = SpectraREWSCalculation(AnSpecParam, TurbsimParam, URef)

    % Analytic spectrum calculation for the wind speed URef
    S       = AnSpecParam.sigma^2 * 4*AnSpecParam.L/URef ./ (1 + 6 * AnSpecParam.f * AnSpecParam.L/URef ).^(5/3);
    
    %Grid Parameters
    kappa               = 12*((AnSpecParam.f/URef).^2+(0.12/AnSpecParam.L).^2).^0.5;
    R                   = 63; %Rotor diameter for the 5MW NREL Turbine
    [Y,Z]               = meshgrid(TurbsimParam.y,TurbsimParam.z-AnSpecParam.h);
    DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
    nPoint              = length(DistanceToHub);
    IsInRotorDisc       = DistanceToHub<=R;
    nPointInRotorDisc   = sum(IsInRotorDisc);

    % loop over every point in the grid
    SUM_gamma_uu       	= zeros(size(AnSpecParam.f));       % allocation
    for iPoint=1:1:nPoint                       % ... all iPoints
        if IsInRotorDisc(iPoint)                %Check if the point is in the area where the rotor is moving or not
            for jPoint=1:1:nPoint               % ... all jPoints
                if IsInRotorDisc(jPoint)        %Check if the point is in the area where the rotor is moving or not
                    Distance        = ((Y(jPoint)-Y(iPoint))^2+(Z(jPoint)-Z(iPoint))^2)^0.5;
                    SUM_gamma_uu    = SUM_gamma_uu + exp(-kappa.*Distance); % coherence is exp(-kappa.*Distance)
                end
            end
         end
    end
    
    % spectra rotor-effective wind speed - rotor model
    S_RR = S/nPointInRotorDisc^2.*SUM_gamma_uu;

return;