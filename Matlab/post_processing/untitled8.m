    %Read the n ROSCO dbg file for URef wind speed and store it in a cell array
    URef = 12;
    DataRosco = ReadRoscoDbgFiles(URef,6);
    DataFAST = ReadFASToutbFiles(URef,6);
    DataRoscoSS = ReadRoscoDbgFiles(URef,6);

%%
%%figure
ax3 = subplot(2,1,1);
hold on; grid on; box on
plot(1:length(DataFAST{1}.BlPitchC1),DataFAST{1}.BlPitchC1)

ax4 = subplot(2,1,2);
hold on; grid on; box on
plot(1:length(DataRosco{1}.WE_b),(180/pi)*DataRosco{1}.WE_b)
plot(1:length(DataRosco{1}.WE_b),(180/pi)*DataRoscoSS{1}.WE_b)
%%
figure 
ax1 = subplot(2,1,1);
hold on; grid on; box on
plot(1:length(DataRosco{1}.GenSpeedF),DataRosco{1}.GenSpeedF)

ax2 = subplot(2,1,2);
hold on; grid on; box on
plot(1:length(DataRosco{1}.RotSpeedF),DataRosco{1}.RotSpeedF)
plot(1:length(DataRosco{1}.RotSpeedF),DataRoscoSS{1}.RotSpeedF)

%%
ax5 = subplot(3,1,2);
hold on; grid on; box on
plot(1:length(DataRosco{1}.WE_Vw),DataRosco{1}.WE_Vw)
plot(1:length(DataRosco{1}.WE_Vw),DataRoscoSS{1}.WE_Vw)

ax6 = subplot(3,1,1);
hold on; grid on; box on
plot(1:length(DataRosco{1}.WE_Vm),DataRosco{1}.WE_Vm)
plot(1:length(DataRosco{1}.WE_Vm),DataRoscoSS{1}.WE_Vm)
%%

ax7 = subplot(2,1,1);
hold on; grid on; box on
plot(1:length(DataRosco{1}.WE_t),DataRosco{1}.WE_t)
plot(1:length(DataRosco{1}.WE_t),DataRoscoSS{1}.WE_t)

%%figure
ax8 = subplot(2,1,2);
hold on; grid on; box on
plot(1:length(DataFAST{1}.GenTq),DataFAST{1}.GenTq)


signal1= detrend(REWS(1,:),'constant');
signal2= detrend(TEREWS(1,:),'constant');
time=0:1/20:600-1/20;
[c,lags]=xcorr(signal1,signal2,'coeff');
figure
hold on; grid on;
plot(lags*1/20,c)
xlim([-1 1]*10)
top = max(c);
xi = interp1(c,lags*1/20,top) ;

figure
hold on; grid on;
plot(time,signal1)
plot(time,signal2)
plot(time+xi,signal2)


    W_cutoff = 0.033333411319;
    W_delay = 2*pi*0.1;
    W_n = W_delay/W_cutoff;
    T_filter = (atan(W_n))/(W_delay);