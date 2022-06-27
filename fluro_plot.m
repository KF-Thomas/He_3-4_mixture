signal = readtable('E:\EXPERIMENT-DATA\NewFile2.csv');
signal = signal.Var1;
background = readtable('E:\EXPERIMENT-DATA\NewFile3.csv');%
background = background.Var1;
%%
t = (0:(size(signal,1)-1)).*2e-6;
stfig('Fluro Sig')
clf
plot(t,signal,'LineWidth',1.9)
hold on
plot(t,background,'LineWidth',1.9)
ylabel('Photodiode Voltage (V)')
xlabel('time (s)')
set(gca,'FontSize',20)
legend('Signal','Background')
%%
stfig('Fluro Sig-B')
clf
plot(t,signal-background)
ylabel('Signal Voltage (V)')
xlabel('time (s)')