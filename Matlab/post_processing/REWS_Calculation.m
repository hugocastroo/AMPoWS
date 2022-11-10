%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a wind field using n different seeds for different wind fields.
%The wind field, has a grid size of 33x33. The REWS and the spectra are calculated and compared with the
%analytic values.
%------------------------------------------------------------
%V1.0 2022.11.09 - Based on Script from David Schlip's lecture
% ----------------------------------
clearvars;clc;%close all;

%% Set the basic info for the simulation
nSeed = 1;                                                                  %Different seeds for different wind fields
Seed_vec = 1:nSeed;
URefMinimum = 4;
URefMaximum = 4;
URefStep = 2;
%%
%Loop for every different wind speed
for URef = URefMinimum:URefStep:URefMaximum
    if(URef < 10)                                                           %%Take just one digit for speeds below 10m/s
        resolution = '%01d';
    else
        resolution = '%02d';
    end
    %Loop for reading every different seed according to the wind speed using readBLGrid function
    velocity = cell(1,nSeed);                                               % allocation
    for iSeed = 1:nSeed
        Seed                = iSeed;
        TurbSimResultFile  	= ['e:\Tesis\Simulationen\Teil1\wind\NTM_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_turbsim.wnd'];
        [velocity{iSeed}, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(TurbSimResultFile);
    end

    %% After reading the files create the analytic spectrum of the wind speed make an estimation for the whole wind field using pwelch.
    % frequency data
    T       = size(velocity{1},1)*dt;
    f_max   = 1/2*1/dt;
    f_min   = 1/T;
    df      = f_min;
    f       = f_min:df:f_max;
    n_f     = length(f);
    
    % SummVars from wnd files
    sigma   = SummVars(4)/100*URef;
    h       = SummVars(1);
    
    % Length scale IEC
    L       = 8.1*42;
    
    % time
    t       = 0:dt:T-dt;
    n_t     = length(t); 
    
    % Analytic spectrum calculation for the wind speed URef
    S       = sigma^2 * 4*L/URef ./ (1 + 6 * f * L/URef ).^(5/3);

    %%Calculation for the whole grid
    kappa               = 12*((f/URef).^2+(0.12/L).^2).^0.5;
    R                   = 63;
    [Y,Z]               = meshgrid(y,z-h);
    DistanceToHub       = (Y(:).^2+Z(:).^2).^0.5;
    nPoint              = length(DistanceToHub);
    IsInRotorDisc       = DistanceToHub<=R;
    nPointInRotorDisc   = sum(IsInRotorDisc);

    % loop over every point in the grid
    SUM_gamma_uu       	= zeros(size(f));       % allocation
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
    
    % estimation
    nBlocks                     = 1;
    SamplingFrequency           = 1/dt;
    n_FFT                       = n_t/nBlocks;
    MyWindow                    = ones(n_FFT,1);
    % MyWindow                    = hamming(n_FFT);
    
    %Create an array with the different spectrum estimations with pwelch using every different wind feld according to the different seeds
    S_est = zeros(nSeed,((n_FFT/2)+1));                                     % allocation                                            
    for iSeed = 1:nSeed
        v_0             = NaN(n_t,1);
        for i_t = 1:1:n_t
            CurrentWind     = squeeze(velocity{iSeed}(i_t,1,:,:));          % extract signal from wind field array and squeze it
            WindField       = CurrentWind(IsInRotorDisc);
            v_0(i_t,1)      = mean(WindField);
        end
        % estimate spectrum
        [S_est(iSeed,:),f_est]   	= pwelch(v_0-mean(v_0),MyWindow,[],n_FFT,SamplingFrequency);
    end

    %Calculate the mean REWS using the array with the different spectrum data from every seed
    S_est_mean          = mean(S_est);

    %plot the result
    figure
    hold on;grid on;box on
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    xlabel('frequency [Hz]','FontSize', 20)
    ylabel('Spectrum [(m/s)^2/Hz]','FontSize', 20)
    plot(f,S_RR,'r','Linewidth',2)                                          %Plot the analytic spectrum
    plot(f_est,S_est_mean,'b','Linewidth',2)                                %Plot the estimation made with the mean values and pwelch

end

%Save plot as pdf for further use if desired
SavePlot = false;                                       
if SavePlot
    pdf = gcf;
    set(pdf,'PaperOrientation','landscape');
    set(pdf,'PaperUnits','normalized');
    set(pdf,'PaperPosition', [0 0 1 1]);
    pdfpath = ['Windspeed_', num2str(URef,resolution), 'ms','_Gridcenter_',num2str(idx_y),'_Seeds',num2str(nSeed)];
    saveas(gcf,pdfpath,'pdf');
end

%% Rosco REWS-Kalman Filter Estimation with the Turbsim REWS
%Plot original
RoscoResultFile  	= ['e:\Tesis\Simulationen\Teil1\sim\1p2_RandSeed1-401_maininput.RO.dbg'];
DataRosco    = ReadROSCOtextIntoStruct(RoscoResultFile);
figure
hold on; grid on; box on
plot(DataRosco.Time,  DataRosco.WE_Vw,'b','Linewidth',2,'DisplayName','Estimated wind speed');
plot(t,v_0,'g','Linewidth',2,'DisplayName','Rotor disc')
xlabel("Time in s");
ylabel("Velocity in m/s");
xlim([0 1200]);
lgd = legend;
%Plot with the signals shifted and the 700 seconds adjustment
figure
hold on; grid on; box on
plot(DataRosco.Time(1:length(DataRosco.Time)-10334)-600,  DataRosco.WE_Vw(1:length(DataRosco.Time)-10334),'b','Linewidth',2,'DisplayName','Estimated wind speed');
plot(DataRosco.Time(5666:end),  DataRosco.WE_Vw(5666:end),'b','Linewidth',2,'DisplayName','Estimated wind speed');
plot(t+4.5,v_0,'g','Linewidth',2,'DisplayName','Rotor disc')
xlabel("Time in s");
ylabel("Velocity in m/s");
xlim([0 1200]);
lgd = legend;

%Compare wind field from 0-600 with 600-1200
%Plot with the signals shifted and the 700 seconds adjustment
figure
hold on; grid on; box on
plot(DataRosco.Time,DataRosco.WE_Vw,'b','Linewidth',2,'DisplayName','Estimated wind speed');
plot(DataRosco.Time-600,  DataRosco.WE_Vw,'g','Linewidth',2,'DisplayName','Estimated wind speed');
xlabel("Time in s");
ylabel("Velocity in m/s");
xlim([0 1200]);
lgd = legend;

% for URef = URefMinimum:URefStep:URefMaximum
%     if(URef < 10)                                                           %%Take just one digit for speeds below 10m/s
%         resolution = '%01d';
%     else
%         resolution = '%02d';
%     end
%     %Loop for reading every different seed according to the wind speed using readBLGrid function
%     DataRosco = cell(1,nSeed);                                               % allocation
%     for iSeed = 1:nSeed
%         Seed                = iSeed;
%         RoscoResultFile  	= ['e:\Tesis\Simulationen\Teil1\sim\Simulationen\1p2_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_maininput.RO.dbg'];
%         DataRosco{iSeed}    = ReadROSCOtextIntoStruct(RoscoResultFile);
%     end
% 
%     ModifiedRoscoWE_Vw = zeros(nSeed,(length(DataRosco{iSeed}.WE_Vw)-1)/20);       % allocation
%     for iSeed = 1:nSeed
%         for i = 1:(length(ModifiedRoscoWE_Vw))
%             ModifiedRoscoWE_Vw(i) = DataRosco{iSeed}.WE_Vw((i*20)-19);
%         end
%     end
%     
%     % frequency data
%     %dt = DataRosco.Time(2)-DataRosco.Time(1);
%     T       = length(ModifiedRoscoWE_Vw)*dt;
%     f_max   = 1/2*1/dt;
%     f_min   = 1/T;
%     df      = f_min;
%     f       = f_min:df:f_max;
%     n_f     = length(f);
%     
%     
%     % time
%     t       = 0:dt:T-dt;
%     n_t     = length(t); 
% 
%     % estimation
%     nBlocks                     = 1;
%     SamplingFrequency           = 1/dt;
%     n_FFT                       = n_t/nBlocks;
%     MyWindow                    = ones(n_FFT,1);
%     % MyWindow                    = hamming(n_FFT);
%     
%     %Create an array with the different spectrum estimations with pwelch using every different wind feld according to the different seeds
%     S_estRosco = zeros(nSeed,((n_FFT/2)+1));                                     % allocation                                            
%     for iSeed = 1:nSeed
%         v_0 = ModifiedRoscoWE_Vw(iSeed,:);
%         % estimate spectrum
%         [S_estRosco(iSeed,:),f_estRosco]   	= pwelch(v_0-mean(v_0),MyWindow,[],n_FFT,SamplingFrequency);
%     end
% 
%     %Calculate the mean REWS using the array with the different spectrum data from every seed
%     S_est_mean          = mean(S_estRosco);
% 
%     plot(f_estRosco,S_est_mean,'g','Linewidth',2)                                %Plot the estimation made with the mean values and pwelch
% end                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   