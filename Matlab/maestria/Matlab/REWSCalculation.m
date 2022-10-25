% -----------------------------
% Script: Reading a bts file and calculating the REWS from the wind field, 
% ------------
% History:

% ----------------------------------
%clear all;clc;close all;
addpath('Windfields')
addpath('dbgFiles')
Data1 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2401_maininput.RO.dbg');
Data2 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2402_maininput.RO.dbg');
Data3 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2403_maininput.RO.dbg');
Data4 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2404_maininput.RO.dbg');
Data5 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2405_maininput.RO.dbg');
Data6 = ReadROSCOtextIntoStruct('1p2_RandSeed1-2406_maininput.RO.dbg');
velocityKF = [Data1.WE_Vw; Data2.WE_Vw; Data3.WE_Vw; Data4.WE_Vw; Data5.WE_Vw; Data6.WE_Vw];
timeKF = [0:0.00625:3600.03125];
clear Data1 Data2 Data3 Data4 Data5 Data6;
% plot(timeKF,velocityKF);
% legend('REWS complete Windfield','TEREWS(KF Turbine)')
%xlim([3050 3600]);

%% Read the wnd files for different URef with different seeds i
addpath('Windfields')
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
for URef = 24:2:24
    velocity = eval(['velocity', int2str(URef),';']);
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
    %URef    = eval(['SummVars', int2str(URef),'(3)';]); % Get URef from the corresponding summvars
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
    tic
    [S_est,f_est]               = pwelch(u-mean(u),hamming(n_FFT),[],n_FFT,SamplingFrequency);
    toc
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
    
    kappa               = 12*((f/URef).^2+(0.12/L).^2).^0.5;
    R                   = 63;
    [Y,Z]               = meshgrid(y,z-h);
    DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
    nPoint              = length(DistanceToHub);
    IsInRotorDisc       = DistanceToHub<=R;
    nPointInRotorDisc   = sum(IsInRotorDisc);
%     tic
%     % loop over ...
%     SUM_gamma_uu       	= zeros(size(f));       % allocation
%     for iPoint=1:1:nPoint                       % ... all iPoints
%         if IsInRotorDisc(iPoint)
%             for jPoint=1:1:nPoint               % ... all jPoints
%                 if IsInRotorDisc(jPoint)
%                     Distance        = ((Y(jPoint)-Y(iPoint))^2+(Z(jPoint)-Z(iPoint))^2)^0.5;
%                     SUM_gamma_uu    = SUM_gamma_uu + exp(-kappa.*Distance);
%                 end
%             end
%          end
%     end
%     toc
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
    tic
    [S_est,f_est]               = pwelch(v_0-mean(v_0),hamming(n_FFT),[],n_FFT,SamplingFrequency);
    toc
    plot(f_est,S_est,'.-')
    plot(f,S_RR,'-')
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    xlabel('frequency [Hz]')
    ylabel('Spectrum [(m/s)^2/Hz]')
    legend(['One point Hub estimate', ['velocity', int2str(URef)]],'Hub analytic','estimate rotor','analytic rotor')
    
    % plot
    figure
    hold on;grid on;box on
    %plot(t,u)
    plot(t,v_0)
    xlabel('time [s]')
    ylabel('wind speed [m/s]')
    %legend('hub height','rotor')
end
%% Compare the analytic spectrum of the whole rotor-effective wind speed with the estimated one using pwelch.
% kappa               = 12*((f/URef).^2+(0.12/L).^2).^0.5;
% R                   = 63;
% [Y,Z]               = meshgrid(y,z-h);
% DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
% nPoint              = length(DistanceToHub);
% IsInRotorDisc       = DistanceToHub<=R;
% nPointInRotorDisc   = sum(IsInRotorDisc);
% 
% % loop over ...
% SUM_gamma_uu       	= zeros(size(f));       % allocation
% for iPoint=1:1:nPoint                       % ... all iPoints
%     if IsInRotorDisc(iPoint)
%         for jPoint=1:1:nPoint               % ... all jPoints
%             if IsInRotorDisc(jPoint)
%                 Distance        = ((Y(jPoint)-Y(iPoint))^2+(Z(jPoint)-Z(iPoint))^2)^0.5;
%                 SUM_gamma_uu    = SUM_gamma_uu + exp(-kappa.*Distance);
%             end
%         end
%      end
% end
% 
% % spectra rotor-effective wind speed
% S_RR = S/nPointInRotorDisc^2.*SUM_gamma_uu;
% 
% % get rotor-effective wind speed
% v_0     = NaN(n_t,1);
% for i_t = 1:1:n_t
%     CurrentWind     = squeeze(velocity(i_t,1,:,:)); 
%     WindField       = CurrentWind(IsInRotorDisc);
%   	v_0(i_t,1)      = mean(WindField);
% end
% 
% % estimation
% nBlocks                     = 32;
% SamplingFrequency           = 1/dt;
% n_FFT                       = n_t/nBlocks;
% [S_est,f_est]               = pwelch(v_0-mean(v_0),hamming(n_FFT),[],n_FFT,SamplingFrequency);
% 
% % plot
% %figure
% % hold on;grid on;box on
% % plot(t,u)
% % plot(t,v_0)
% % xlabel('time [s]')
% % ylabel('wind speed [m/s]')
% % legend('hub height','rotor')
% 
% % plot
% %figure
% hold on;grid on;box on
% plot(f_est,S_est,'.-')
% plot(f,S_RR,'-')
% set(gca,'xScale','log')
% set(gca,'yScale','log')
% xlabel('frequency [Hz]')
% ylabel('Spectrum [(m/s)^2/Hz]')
% legend('estimate','analytic')

%% e) Compare the analytic coherence with an estimated one from two points at hub height with a distance of 20m using mscohere.

% get signals
idx_y_1             = 17;
idx_y_2             = 22;
idx_z_1             = 17;
idx_z_2             = 17;
u_1                 = velocity(:,1,idx_y_1,idx_z_1);
u_2                 = velocity(:,1,idx_y_2,idx_z_2);
Distance            = ((y(idx_y_1)-y(idx_y_2))^2+(z(idx_z_1)-z(idx_z_2))^2)^0.5;
fprintf('Distance: %4.1f m.\n',Distance)

% coherence gamma = exp(-kappa.*Distance)
gamma               = exp(-kappa.*Distance);

% estimate coherence
nBlocks           	= 16;
SamplingFrequency   = 1/dt;
n_FFT           	= n_t/nBlocks;
[gamma_Sq_est,f_est]= mscohere(u_1-mean(u_1),u_2-mean(u_2),hamming(n_FFT),[],n_FFT,SamplingFrequency);

figure
hold on;grid on;box on
plot(f_est,gamma_Sq_est,'-')
plot(f,gamma.^2,'-')
set(gca,'xScale','log')
xlabel('frequency [Hz]')
ylabel('Squared Coherence [-]')
legend('estimate','analytic')
