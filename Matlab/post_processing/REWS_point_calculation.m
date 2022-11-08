%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a single point in a wind field from n different simulations, and n different seeds.
%In this case, the point in the wind field, is the point in the center of
%wind field (33x33), the point in front o the hub. The calculation
%ilustrates the wind speed at Hub height since the point is in front of the
%hub. The REWS and the spectra are calculated and compared with the
%analytic values. After the calculation of a single point, the REWS of the
%whole grid is calculated.
%------------------------------------------------------------
%V1.0 2022.11.07 - Based on Script from David Schlip's lecture
% ----------------------------------
clearvars;clc;%close all;

%% Read the wind field into Matlab with readBLgrid.m.
%Different seeds for different wind fields
nSeed = 6;      
Seed_vec = 1:nSeed;
%Loop for the wind speeds
for URef = 4:2:24
    %%Take just one digit for speeds below 10m/s
    if(URef < 10)
        resolution = '%01d';
    else
        resolution = '%02d';
    end
    %Loop for reading every different seed according to the wind speed
    for iSeed = 1:nSeed
        Seed                = iSeed;
        TurbSimResultFile  	= ['e:\Tesis\Simulationen\Teil1\wind\NTM_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_turbsim.wnd'];
        [velocity{iSeed}, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(TurbSimResultFile);
    end

    %% After reading the files create the analytic spectrum of a single point and an estimated one for the wind at hub height using pwelch.
    
    % frequency data
    T       = size(velocity{1},1)*dt;
    f_max   = 1/2*1/dt;
    f_min   = 1/T;
    df      = f_min;
    f       = [f_min:df:f_max];
    n_f     = length(f);
    
    % SummVars from wnd files
    sigma   = SummVars(4)/100*URef;
    
    % Length scale IEC
    L       = 8.1*42;
    
    % time
    t       = 0:dt:T-dt;
    n_t     = length(t); 
    
    % Analytic spectrum calculation
    S       = sigma^2 * 4*L/URef ./ (1 + 6 * f * L/URef ).^(5/3);
    
    % modified analytic spectrum to have the right std deviation
    S_mod   = S/(sum(S)*df)*sigma^2;
    
    %% estimation calculation using just one block, since it will be meaned over the different seeds
    n_Blocks                 	= 1;
    SamplingFrequency           = 1/dt;
    n_FFT                       = n_t/n_Blocks;
    MyWindow                    = ones(n_FFT,1);
    % MyWindow                    = hamming(n_FFT);
    
    % idx for grid point - Hub Height point
    idx_y 	= 17;
    idx_z  	= 17;
    
    %Create an array with the different spectrum estimations with pwelch using
    %the every different wind feld according to the different seeds, JUST FOR
    %ONE POINT OF THE GRID
    for iSeed = 1:nSeed
        % extract signal
        u                           = velocity{iSeed}(:,1,idx_y,idx_z);
    
        % estimate spectrum
        [S_est(iSeed,:),f_est]   	= pwelch(u-mean(u),MyWindow,[],n_FFT,SamplingFrequency);
    end
    
    % calculate mean spectrum
    S_est_mean          = mean(S_est);
    
    % plot analytic and mean
    figure
    hold on;grid on;box on
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    xlabel('frequency [Hz]','FontSize', 20)
    ylabel('Spectrum [(m/s)^2/Hz]','FontSize', 20)
    %ylim([1e-2 1e3])
    %xlim([1e-3 1e1])
    plot(f,S,'r','Linewidth',1.5)
    plot(f,S_mod,'b','Linewidth',1.5)
    plot(f_est,S_est_mean,'k','Linewidth',1.5)
    
    % plot each single signal if desired
    PlotSingleSpectra           = false;
    
    % plot single spectrum
    if PlotSingleSpectra
        for iSeed = 1:nSeed
            plot(f_est,S_est(iSeed,:),'.-')
        end
        LegendCell = cellstr(strcat('estimateSeed',num2str(Seed_vec','%03d')));
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
end