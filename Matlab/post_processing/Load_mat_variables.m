%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: Load the the variables from the mat file to use them in the
%simulation
%------------------------------------------------------------
%V1.0 2022.10.28
%------------------------------------------------------------
function [velocity, y, z, nz, ny, dz, dy, dt, z1, SummVars] = Load_mat_variables(WindSpeed,i)
    for URef = WindSpeed
        for i = i %Loop for the seeds
            resolution = '%02d';
            if(URef < 10)
                resolution = '%01d';
            end
            load(['NTM_RandSeed1-',num2str(URef,resolution),num2str(i,'%02d'),'_turbsim.mat'],'velocity', 'y', 'z', 'nz', 'ny', 'dz', 'dy', 'dt', 'z1', 'SummVars');
        end
    end
end