%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script uses the function readBLgrid to read n wind fields
%generated in turbsim and return the velocities of the wind fields in a
%cell array and the information data about the files in a struct package
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------

function [velocity, TurbSimParam] = ReadWndFiles(URef,nSeed)
%

    if(URef < 10)                                                           %%Take just one digit for speeds below 10m/s
        resolution = '%01d';
    else
        resolution = '%02d';
    end

    velocity = cell(1,nSeed);                                               % allocation
    for iSeed = 1:nSeed
        Seed                = iSeed;
        TurbSimResultFile  	= ['e:\Tesis\Simulationen\Teil1\wind\NTM_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_turbsim.wnd'];
        [velocity{iSeed}, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(TurbSimResultFile);
        TurbSimParam.y = y;
        TurbSimParam.z = z;
        TurbSimParam.nz = nz;
        TurbSimParam.ny = ny;
        TurbSimParam.dz = dz;
        TurbSimParam.dy = dy;
        TurbSimParam.dt = dt;
        TurbSimParam.zHub = zHub;
        TurbSimParam.z1 = z1;
        TurbSimParam.SummVars = SummVars;
    end

return;