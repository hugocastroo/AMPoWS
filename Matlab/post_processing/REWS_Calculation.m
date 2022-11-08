%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script calculates the REWS of a single point in a wind field.
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

%% b) Read the wind field into Matlab with readBLgrid.m.
nSeed = 6;
Seed_vec = 1:nSeed;
for URef = 4:2:4 %Loop for the wind speeds
    if(URef < 10)
        resolution = '%01d';
    else
        resolution = '%02d';
    end

    for iSeed = 1:nSeed
        Seed                = iSeed;
        TurbSimResultFile  	= ['e:\Tesis\Simulationen\Teil1\wind\NTM_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_turbsim.wnd'];
        [velocity{iSeed}, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(TurbSimResultFile);
    end

    %% Compare the analytic spectrum of a single point with an estimated one for the wind at hub height using pwelch.
    % frequency
    T       = size(velocity{1},1)*dt;
    f_max   = 1/2*1/dt;
    f_min   = 1/T;
    df      = f_min;
    f       = [f_min:df:f_max];
    n_f     = length(f);
    
    % SummVars
    sigma   = SummVars(4)/100*URef;
    h       = SummVars(1);
    
    % Length scale IEC
    L       = 8.1*42;
    
    % time
    t       = 0:dt:T-dt;
    n_t     = length(t); 
    
    % spectrum
    S       = sigma^2 * 4*L/URef ./ (1 + 6 * f * L/URef ).^(5/3);

    %%Calculation for the whole grid
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
        CurrentWind     = squeeze(velocity{iSeed}(i_t,1,:,:)); 
        WindField       = CurrentWind(IsInRotorDisc);
        v_0(i_t,1)      = mean(WindField);
    end
    
    % estimation
    nBlocks                     = 1;
    SamplingFrequency           = 1/dt;
    n_FFT                       = n_t/nBlocks;
    MyWindow                    = ones(n_FFT,1);
    % MyWindow                    = hamming(n_FFT);
    
    % extract signal from wind field array and squeze it
    for iSeed = 1:nSeed
        v_0     = NaN(n_t,1);
        for i_t = 1:1:n_t
            CurrentWind     = squeeze(velocity{iSeed}(i_t,1,:,:)); 
            WindField       = CurrentWind(IsInRotorDisc);
            v_0(i_t,1)      = mean(WindField);
        end
    
        % estimate spectrum
        [S_est(iSeed,:),f_est]   	= pwelch(v_0-mean(v_0),MyWindow,[],n_FFT,SamplingFrequency);
    end
    
    %Calculate the mean REWS using the array with the different spectrum
    %data from every seed
    S_est_mean          = mean(S_est);

    %plot the result
    figure
    hold on;grid on;box on
    set(gca,'xScale','log')
    set(gca,'yScale','log')
    xlabel('frequency [Hz]','FontSize', 20)
    ylabel('Spectrum [(m/s)^2/Hz]','FontSize', 20)
    plot(f,S_RR,'m','Linewidth',2)
    plot(f_est,S_est_mean,'y','Linewidth',2)

end

SavePlot = false;
if SavePlot
    pdf = gcf;
    set(pdf,'PaperOrientation','landscape');
    set(pdf,'PaperUnits','normalized');
    set(pdf,'PaperPosition', [0 0 1 1]);
    pdfpath = ['Windspeed_', num2str(URef,resolution), 'ms','_Gridcenter_',num2str(idx_y),'_Seeds',num2str(nSeed)];
    saveas(gcf,pdfpath,'pdf');
end