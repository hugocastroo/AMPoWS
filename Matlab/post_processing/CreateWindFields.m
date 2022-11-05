%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: Create wnd files for the wind fields if needed for n seeds
%------------------------------------------------------------
%V1.0 2022.10.28
%------------------------------------------------------------
function CreateWindFields(windSpeed,seeds)
    for URef = windSpeed %Loop for the wind speeds
        resolution = '%02d';
        if(URef < 10)
            resolution = '%01d';
        end
        for i = 1:1:seeds %Loop for the seeds
            OutputFile  = ['NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim.wnd'];
            cd WindFiles\
            if ~exist(OutputFile,'file')
                dos(['TurbSim_x64.exe ','NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim.inp' ]);
            end
            cd ..
        end
    end
end