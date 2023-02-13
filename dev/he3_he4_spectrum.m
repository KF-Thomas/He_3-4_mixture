temp = readtable('NewFile3.csv');
stfig('SFP Spectrum for mixed light');
clf
y_vals = smooth(temp{:,2},8);
y_vals = y_vals-mean(y_vals(1:10));
plot(y_vals./max(y_vals),'linewidth',2)
ylabel('Amplitude (arb. u.)')
xlabel('Frequency (arb. u.)')
set(gca,'FontSize',17)
set(gca,'XTickLabel',[]);

ylim([0 1.05])
