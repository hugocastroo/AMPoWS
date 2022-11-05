%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: Loads the wnd turbsim files and store the variables in mat files
%------------------------------------------------------------
%V1.0 2022.10.28
%------------------------------------------------------------
function Load_wnd_files(windSpeed,seeds)
    for URef = windSpeed %Loop for the wind speeds
        resolution = '%02d';
        if(URef < 10)
            resolution = '%01d';
        end
        for i = 1:1:seeds %Loop for the seeds
            if(~exist(['NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim.mat'],'file'))
                %Read the variables from the wnd file
                [velocity, y, z, nz, ny, dz, dy, dt, ~, z1, SummVars] = readBLgrid(['NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim.wnd']);
                %Store the variables in a mat file to be able to use them later
                save(['NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim'],'velocity', 'y', 'z', 'nz', 'ny', 'dz', 'dy', 'dt', 'z1', 'SummVars');
            end
        end
    end
end