% -----------------------------
% Script: REWS comparison - ROSCO KF and TurbSim rotor field
% ------------
%%clear all;clc;close all;
addpath('Windfields')
addpath('dbgFiles')
%% Read the data from the KF dbg files
for URef = 4:2:24
    for i = 1:6
        resolution = '%02d';
        if(URef < 10)
            resolution = '%01d';
        end
        Data = ReadROSCOtextIntoStruct(['1p2_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_maininput.RO.dbg']);
        if(i == 1)
            eval(['velocityKF', int2str(URef),' = Data.WE_Vw(1:end-1);']);
        else
            velocityKFAux = [eval(['velocityKF', int2str(URef)]); Data.WE_Vw(1:end-1)];
            eval(['velocityKF', int2str(URef),' = velocityKFAux;']);
        end
    end
end
timeKF = [0:0.00625:3600-0.00625];
clear velocityKFAux Data
%Data1 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2401_maininput.RO.dbg');
%Concatenate the KF velocities
%velocityKF = [Data1.WE_Vw(1:end-1); Data2.WE_Vw(1:end-1); Data3.WE_Vw(1:end-1); Data4.WE_Vw(1:end-1); Data5.WE_Vw(1:end-1); Data6.WE_Vw(1:end-1)];
%timeKF = [0:0.00625:3600-0.00625];
%Clean the variables for memory
%clear Data1 Data2 Data3 Data4 Data5 Data6;

%% Read the wnd files for different URef with different seeds i
for URef = 4:2:24
    for i = 1:6
        resolution = '%02d';
        if(URef < 10)
            resolution = '%01d';
        end
        [velocityMain, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(['NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim.wnd']);
        if(i == 1)
            eval(['velocity', int2str(URef),' = velocityMain;']);
            eval(['SummVars', int2str(URef),' = SummVars;']);
        else
            velocityAux = [eval(['velocity', int2str(URef)]); velocityMain];
            eval(['velocity', int2str(URef),' = velocityAux;']);
        end
    end
end
%% Calculate the analytic spectrum of a single point at hub height and compare it with an stimated using pwelch
for URef = 4:2:24
    velocity = eval(['velocity', int2str(URef),';']);
    velocityKF = eval(['velocityKF', int2str(URef),';']);
    % extract signal - point in the coordinates 17,17, Point that is exactly in
    % the rotor height JUST ONE POINT
    idx_y               = 17;
    idx_z               = 17;
    u                   = velocity(:,1,idx_y,idx_z);
    % frequency data
    T       = size(velocity,1)*dt;
    f_max   = 1/2*1/dt;
    f_min   = 1/T;
    df      = f_min;
    f       = [f_min:df:f_max];
    n_f     = length(f);
    
    % SummVars Get information from the file variables
    sigma   = eval(['SummVars', int2str(URef),'(4)';])/100*URef;
    h       = eval(['SummVars', int2str(URef),'(1)';]);
    
    % Length scale IEC
    L       = 8.1*42;
    
    % time
    t       = 0:dt:T-dt;
    n_t     = length(t); 
    
    % spectrum calculation
    S       = sigma^2 * 4*L/URef ./ (1 + 6 * f * L/URef ).^(5/3);
    
    % estimation "More blocks mean better resolution, but less frequency lenght since the time series is beeing splitted into more parts"
    nBlocks                     = 32;
    SamplingFrequency           = 1/dt;
    n_FFT                       = n_t/nBlocks;
    [S_est,f_est]               = pwelch(u-mean(u),hamming(n_FFT),[],n_FFT,SamplingFrequency);

    % plot
    figure
    hold on;grid on;box on
    plot(f_est,S_est,'.-')
    plot(f,S,'-')
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    xlabel('frequency [Hz]')
    ylabel('Spectrum [(m/s)^2/Hz]')
    legend(['One point Hub estimate', ['velocity', int2str(URef)]],'Hub analytic')
    %--------------------Finishes one point estimate
    kappa               = 12*((f/URef).^2+(0.12/L).^2).^0.5;
    R                   = 63;
    [Y,Z]               = meshgrid(y,z-h);
    DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
    nPoint              = length(DistanceToHub);
    IsInRotorDisc       = DistanceToHub<=R;
    nPointInRotorDisc   = sum(IsInRotorDisc);
    % loop over ...
   SUM_gamma_uu       	= zeros(size(f));       % allocation
    for iPoint=1:1:nPoint                       % ... all iPoints
        if IsInRotorDisc(iPoint)
            for jPoint=1:1:nPoint               % ... all jPoints
                if IsInRotorDisc(jPoint)
                    Distance        = ((Y(jPoint)-Y(iPoint))^2+(Z(jPoint)-Z(iPoint))^2)^0.5;
                    SUM_gamma_uu    = SUM_gamma_uu + exp(-kappa.*Distance);
                end
            end
         end
    end
    % spectra rotor-effective wind speed
    
    S_RR = S/nPointInRotorDisc^2.*SUM_gamma_uu;
    
    % get rotor-effective wind speed
    v_0     = NaN(n_t,1);
    for i_t = 1:1:n_t
        CurrentWind     = squeeze(velocity(i_t,1,:,:)); 
        WindField       = CurrentWind(IsInRotorDisc);
  	    v_0(i_t,1)      = mean(WindField);
    end
    
    % estimation
    nBlocks                     = 32;
    SamplingFrequency           = 1/dt;
    n_FFT                       = n_t/nBlocks;
    [S_est,f_est]               = pwelch(v_0-mean(v_0),hamming(n_FFT),[],n_FFT,SamplingFrequency);

    %PSD KF
    nBlocks                     = 32;
    SamplingFrequency           = 1/(timeKF(2)-timeKF(1));
    n_t     = length(timeKF);
    n_FFT                       = n_t/nBlocks;
    [S_estKF,f_estKF]               = pwelch(velocityKF-mean(velocityKF),hamming(n_FFT),[],n_FFT,SamplingFrequency);
    
    %PSD plots Estimate Rotor, Rotor Analytic, Estimate KF
    plot(f_est,S_est,'.-')
    plot(f,S_RR,'-')
    plot(f_estKF,S_estKF,'.-')
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    legend(['One point Hub estimate', ['velocity', int2str(URef)]],'Hub analytic','estimate rotor','analytic rotor','TEREWS(KF)')
    xlabel('frequency [Hz]')
    ylabel('Spectrum [(m/s)^2/Hz]')

    %Plots Wind Fields
    figure
    hold on;grid on;box on
    %%plot(t,u)
    plot(t,v_0)
    plot(timeKF,velocityKF);
    xlabel('time [s]')
    ylabel('wind speed [m/s]')
    legend(['REWS complete Windfield', int2str(URef)],'TEREWS(KF)')
end
%% Analytic TEREWS - At the end is the same as the analytic REWS BUT IT HAS A LONGER SPECTRUM SINCE THE ROSCO CONTROLLER IS SAMPLING AT HIGHER FREQUENCIES
% frequency data

dtKF = timeKF(2)-timeKF(1);
TKF       = length(velocityKF)*dtKF;
f_maxKF   = 1/2*1/dtKF;
f_minKF   = 1/TKF;
dfKF      = f_minKF;
fKF       = [f_minKF:dfKF:f_maxKF];
SKF       = sigma^2 * 4*L/URef ./ (1 + 6 * fKF * L/URef ).^(5/3);
kappaKF               = 12*((fKF/URef).^2+(0.12/L).^2).^0.5;
R                   = 63;
[Y,Z]               = meshgrid(y,z-h);
DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
nPoint              = length(DistanceToHub);
IsInRotorDisc       = DistanceToHub<=R;
nPointInRotorDisc   = sum(IsInRotorDisc);
% loop over ...
SUM_gamma_KF       	= zeros(size(fKF));       % allocation
for iPoint=1:1:nPoint                       % ... all iPoints
    if IsInRotorDisc(iPoint)
        for jPoint=1:1:nPoint               % ... all jPoints
            if IsInRotorDisc(jPoint)
                Distance        = ((Y(jPoint)-Y(iPoint))^2+(Z(jPoint)-Z(iPoint))^2)^0.5;
                SUM_gamma_KF    = SUM_gamma_KF + exp(-kappaKF.*Distance);
            end
        end
     end
end
% spectra rotor-effective wind speed

S_RRKF = SKF/nPointInRotorDisc^2.*SUM_gamma_KF;

figure
%PSD plots Estimate Rotor, Rotor Analytic, Estimate KF
hold on;grid on;box on
plot(f_estKF,S_estKF,'.-')
plot(fKF,S_RRKF,'-')
plot(f_est,S_est,'.-')
plot(f,S_RR,'-')
set(gca,'xScale','log')
set(gca,'yScale','log')
legend('TEREWS(KF)','Analytic TEREWS(KF)','estimate rotor','analytic rotor')
xlabel('frequency [Hz]')
ylabel('Spectrum [(m/s)^2/Hz]')
%% Coherence

%% e) Compare the analytic coherence with an estimated one from two points at hub height with a distance of 20m using mscohere.

dtKF = timeKF(2)-timeKF(1);
TKF       = length(velocityKF)*dtKF;
f_maxKF   = 1/2*1/dtKF;
f_minKF   = 1/TKF;
dfKF      = f_minKF;
fKF       = [f_minKF:dfKF:f_maxKF];
SKF       = sigma^2 * 4*L/URef ./ (1 + 6 * fKF * L/URef ).^(5/3);
kappaX               = 12*((fKF/URef).^2+(0.12/L).^2).^0.5;
SUM_gamma_X       	= zeros(size(fKF));       % allocation

    for jPoint=1:1:nPoint               % ... all jPoints
        if IsInRotorDisc(jPoint)
            DistanceX        = ((Y(jPoint))^2+(Z(jPoint))^2)^0.5;
            SUM_gamma_X    = SUM_gamma_X + exp(-kappaX.*DistanceX);
        end
    end
S_X = SKF/nPointInRotorDisc.*SUM_gamma_X;
cohsq= zeros(1,8192);       % allocation
for i = 1:8192
    cohsq(i) = abs(S_X(35*i))^2/(S_RRKF(35*i)*S_RR(i));
end
grl=(zeros(1,length(S_X))); % allocation
for i = 1:length(S_X)
    grl(i) = abs((S_X)/(S_RRKF));
end

T       = length(cohsq)*dt;
f_max   = 1/dt;
f_min   = 1/T;
df      = f_min;
fcoh    = [f_min:df:f_max];

figure
hold on;grid on;box on
plot(fKF,((SUM_gamma_X/nPointInRotorDisc).^2))
plot(f,((SUM_gamma_uu/nPointInRotorDisc^2).^2))
plot(fKF,((SUM_gamma_KF/nPointInRotorDisc^2).^2))
plot(fcoh,cohsq)
set(gca,'xScale','log')
xlabel('frequency [Hz]')
ylabel('Squared Coherence [-]')
legend('CrossSpectra','RotorSpectra','KFSpectra','CoherenceByHand')

figure
hold on;grid on;box on
plot(fKF,grl)
set(gca,'xScale','log')
xlabel('frequency [Hz]')
ylabel('Grl [-]')
%% Zeitversatz according to the 0.20944 F_WECornerFreq - Corner frequency (-3dB point) in the first order low pass filter for the wind speed estimate [rad/s].

    %Plots Wind Fields
    figure
    hold on;grid on;box on
    %%plot(t,u)
    plot(t,v_0)
    plot(timeKF,velocityKF);
    plot((timeKF-(1/(2*pi*0.20944))),velocityKF);
    xlabel('time [s]')
    ylabel('wind speed [m/s]')
    title("REWS-TEREWS Time offset 1/(2*pi*0.20944)");
    legend(['REWS Windfeld', int2str(URef)],'TEREWS(KF) Original','TEREWS(KF) Offset')
    xlim([3500 3600]);