%% Post Processing OpenFAST V2.5 and V3.2.1 Simulations with ROSCO Controller - Master Thesis Hugo Valentin Castro Saenz
%------------------------------------------------------------
%Script: This script uses the function ReadROSCOtextIntoStruct to read n
%debug files from the ROSCO controller - OpenFAST simulation The n
%structures are stored in a cell array
%------------------------------------------------------------
%V1.0 2022.11.12 - HC
% ----------------------------------

function DataRosco = ReadRoscoDbgFiles(URef,nSeed)

    if(URef < 10)                                                           %%Take just one digit for speeds below 10m/s
        resolution = '%01d';
    else
        resolution = '%02d';
    end

    DataRosco = cell(1,nSeed);                                               % allocation
    for iSeed = 1:nSeed
        Seed                = iSeed;
        if(URef <= 16)                                                           %%Take just one digit for speeds below 10m/s
            RoscoResultFile  	= ['e:\Tesis\Simulationen\Teil1\sim\1p2_',num2str(URef,resolution),'_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_maininput.RO.dbg'];
        else
            RoscoResultFile  	= ['e:\Tesis\Simulationen\Teil1\sim\1p_',num2str(URef,resolution),'_RandSeed1-',num2str(URef,resolution),num2str(Seed,'%02d'),'_maininput.RO.dbg'];
        end
        
        DataRosco{Seed}    = ReadROSCOtextIntoStruct(RoscoResultFile);
    end

return;