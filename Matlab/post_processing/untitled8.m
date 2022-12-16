%     %Read the n ROSCO dbg file for URef wind speed and store it in a cell array
%     URef = 12;
%     DataRosco = ReadRoscoDbgFiles(URef,6);
%     DataFAST = ReadFASToutbFiles(URef,6);
%     DataRoscoSS = ReadRoscoDbgFiles(URef,6);

RoscoResultFile1  	= 'c:\OpenFAST\openfast\reg_tests\r-test\glue-codes\openfast\5MW_Land_DLL_WTurb\5MW_Land_DLL_WTurb.outb';

DataFAST = cell(1,2);
DataFAST{1}    = ReadFASTbinaryIntoStruct(RoscoResultFile1);
%%
%figure
ax1 = subplot(6,1,1);
hold on; grid on; box on
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.BldPitch1,DisplayName='LossTorqueWithAdjustment')
%plot(DataFAST{1, 2}.Time,DataFAST{1, 2}.BldPitch1)
title('Pitch')
legend;

ax2 = subplot(6,1,2);
hold on; grid on; box on
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.RotSpeed)
%plot(DataFAST{1, 2}.Time,DataFAST{1, 2}.RotSpeed)
title('RotSpeed')

ax3 = subplot(6,1,3);
hold on; grid on; box on
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.GenSpeed)
%plot(DataFAST{1, 2}.Time,DataFAST{1, 2}.GenSpeed)
title('GenSpeed')

ax4 = subplot(6,1,4);
hold on; grid on; box on
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.RotTorq*1000)
%plot(DataFAST{1, 2}.Time,DataFAST{1, 1}.LSSTipMzs)
title('RotTrq')

ax5 = subplot(6,1,5);
hold on; grid on; box on
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.GenTq)
%plot(DataFAST{1, 2}.Time,DataFAST{1, 2}.GenTq)
%yline(43.093,'HandleVisibility','off');
title('GenTrq')

ax6 = subplot(6,1,6);
hold on; grid on; box on
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.GenPwr)
%plot(DataFAST{1, 2}.Time,DataFAST{1, 2}.GenPwr)
title('GenPwr')
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x');


%%
% figure 
% ax1 = subplot(2,1,1);
% hold on; grid on; box on
% plot(1:length(DataRosco{1}.GenSpeedF),DataRosco{1}.GenSpeedF)
% 
% ax2 = subplot(2,1,2);
% hold on; grid on; box on
% plot(1:length(DataRosco{1}.RotSpeedF),DataRosco{1}.RotSpeedF)
% plot(1:length(DataRosco{1}.RotSpeedF),DataRoscoSS{1}.RotSpeedF)
% 
% %%
% ax5 = subplot(3,1,2);
% hold on; grid on; box on
% plot(1:length(DataRosco{1}.WE_Vw),DataRosco{1}.WE_Vw)
% plot(1:length(DataRosco{1}.WE_Vw),DataRoscoSS{1}.WE_Vw)
% 
% ax6 = subplot(3,1,1);
% hold on; grid on; box on
% plot(1:length(DataRosco{1}.WE_Vm),DataRosco{1}.WE_Vm)
% plot(1:length(DataRosco{1}.WE_Vm),DataRoscoSS{1}.WE_Vm)
% %%
% 
% ax7 = subplot(2,1,1);
% hold on; grid on; box on
% plot(1:length(DataRosco{1}.WE_t),DataRosco{1}.WE_t)
% plot(1:length(DataRosco{1}.WE_t),DataRoscoSS{1}.WE_t)
% 
% %%figure
% ax8 = subplot(2,1,2);
% hold on; grid on; box on
% plot(1:length(DataFAST{1}.GenTq),DataFAST{1}.GenTq)
% 
% 
% signal1= detrend(REWS(1,:),'constant');
% signal2= detrend(TEREWS(1,:),'constant');
% time=0:1/20:600-1/20;
% [c,lags]=xcorr(signal1,signal2,'coeff');
% figure
% hold on; grid on;
% plot(lags*1/20,c)
% xlim([-1 1]*10)
% top = max(c);
% xi = interp1(c,lags*1/20,top) ;
% 
% figure
% hold on; grid on;
% plot(time,signal1)
% plot(time,signal2)
% plot(time+xi,signal2)
% 
% 
%     W_cutoff = 0.033333411319;
%     W_delay = 2*pi*0.1;
%     W_n = W_delay/W_cutoff;
%     T_filter = (atan(W_n))/(W_delay);

figure
hold on; grid on; box on

plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTDxt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTDyt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTDzt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTDxi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTDyi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTDzi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTVxt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTVyt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTVzt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTVxi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTVyi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTVzi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTAxt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTAyt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTAzt)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTAxi)
plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTAyi)

    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmTAzi)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRDxi)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRDyi)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRDzi)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRVxt)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRVyt)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRVzt)
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRVxi)
    
    plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRVyi)

plot(DataFAST{1, 1}.Time,(DataFAST{1, 1}.PtfmTVzi))

%%
% figure 
% hold on; grid on
% plot(DataFAST{1, 1}.Time,DataFAST{1, 1}.PtfmRVyi)
% yline(max(DataFAST{1, 1}.PtfmRVyi),'--',num2str(max(DataFAST{1, 1}.PtfmRVyi)),'HandleVisibility','on');
% yline(min(DataFAST{1, 1}.PtfmRVyi),'--',num2str(min(DataFAST{1, 1}.PtfmRVyi)),'HandleVisibility','on');
% yline(mean(DataFAST{1, 1}.PtfmRVyi),'--',num2str(mean(DataFAST{1, 1}.PtfmRVyi)),'HandleVisibility','on');