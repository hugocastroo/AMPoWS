%Steady states calculations for different wind speeds
clc
SteadyStates = load('SteadyStatesNREL5MW_FBNREL_SLOW','v_0','Omega','theta','x_T','M_g');

for URef=4:2:24

    theta          	= interp1(SteadyStates.v_0,SteadyStates.theta   ,URef,'linear','extrap');
    Omega          	= interp1(SteadyStates.v_0,SteadyStates.Omega   ,URef,'linear','extrap');
    x_T                = interp1(SteadyStates.v_0,SteadyStates.x_T     ,URef,'linear','extrap');
    
    %fprintf(['Wind Speed: ',num2str(URef),'\n'])
    %fprintf(['BlPitch: ', num2str((180/pi)*theta),'\n'])
    %fprintf(['RotSpeed: ', num2str(radPs2rpm(Omega)),'\n'])
    fprintf(['', num2str(x_T),'\n'])

end
