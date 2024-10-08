% plot LTE-A comparison curves

clear all
close all
clc

% open LTE-A reference figure

cd('../results');

% load parameters
nCQI = 15;
load('5GSimReferenceCurves');
load('LTEASimReferenceCurves');
tempLegend = cell(nCQI+2,1);

% throughput plot
figure();
linecolors = colormap('lines');
for iCQI = 1:nCQI
    plot(throughputResult{iCQI,2},throughputResult{iCQI,1}*1e-3, 'color',linecolors(iCQI,:));
    hold on;
    tmpPlot = plot(squeeze(throughputOverall(iCQI,1,:)),squeeze(throughputOverall(iCQI,2,:))/1e6, 'color',linecolors(iCQI,:), 'marker', 'square');
    set(get(get(tmpPlot,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    tempLegend{iCQI} = ['CQI',num2str(iCQI)];
end
grid on;
plot(0,0,'color','black','linestyle','-','Visible','on');
plot(0,0,'color','black','linestyle','-','marker','square','Visible','on');
tempLegend{nCQI+1} = 'LTE-A Sim';
tempLegend{nCQI+2} = '5G Sim';
xlabel('SNR in dB');
ylabel('Throughput in MBit/s');
xlim([-10 35]);
legend(tempLegend,'location','southeast');

% print figure
pause(1);
set(gcf, 'units', 'centimeters');
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6 4]);
print('-dpdf','throughput_comparison');

% FER plot
figure();
linecolors = colormap('lines');
for iCQI = 1:nCQI
    plot(BLERResult{iCQI,2},BLERResult{iCQI,1}, 'color',linecolors(iCQI,:));
    hold on;
    tmpPlot = plot(squeeze(FEROverall(iCQI,1,:)),squeeze(FEROverall(iCQI,2,:)), 'color',linecolors(iCQI,:), 'marker', 'square');
	set(get(get(tmpPlot,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
h = get(gcf,'CurrentAxes');
set(h,'YScale','log');
grid on;
plot(0,0,'color','black','linestyle','-','Visible','on');
plot(0,0,'color','black','linestyle','-','marker','square','Visible','on');
xlabel('SNR in dB');
ylabel('Frame Error Ratio');
xlim([-10 35]);
legend(tempLegend,'location','southeast');

% print figure
pause(1);
set(gcf, 'units', 'centimeters');
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 6 4]);
print('-dpdf','FER_comparison');