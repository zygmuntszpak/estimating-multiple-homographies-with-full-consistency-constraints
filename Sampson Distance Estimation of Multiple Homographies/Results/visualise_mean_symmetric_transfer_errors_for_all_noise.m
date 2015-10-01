%   Function: visualise_mean_symmetric_transfer_errors_for_all_noise
%
%   Plots lines detailing the overall performance of various homography 
%   estimation methods for various noise levels.  
%
%   Parameters:
%
%      listOfMeanErrorsForEachMethodAndSimulation 	- a struct containing the mean root mean square
%													  error of various homography estimation methods
%													  for different noise levels.
%
%   Returns: 
% 
%
%   See Also:
%
%  compute_symmetric_transfer_error_for_all_methods
%  visualise_symmetric_transfer_errors_for_fixed_noise
%  
%
%  Zygmunt L. Szpak (c) 2012
%  Last modified 15/5/2012 

function visualise_mean_symmetric_transfer_errors_for_all_noise( listOfMeanErrorsForEachMethodAndSimulation )

numberOfNoiseLevels = length(listOfMeanErrorsForEachMethodAndSimulation);

myFig = figure;

hold on
% starting noise level
% this should probably not be hard-coded here, but passed into the function
noiseLevel = 0.2;

yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% NDLT
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{NDLT_CONST}]
end
plot(yData,xData,'Color',[0.752941176470588 0.313725490196078 0.301960784313725],'LineWidth',2,'LineStyle',':','MarkerSize',6,'Marker','o','MarkerEdgeColor', [0.752941176470588 0.313725490196078 0.301960784313725])


yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% FNS
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{FNS_CONST}]
end
plot(yData,xData,'Color',[0.607843137254902  0.733333333333333 0.349019607843137],'LineWidth',2,'LineStyle','-.','MarkerSize',6,'Marker','s','MarkerEdgeColor', [0.607843137254902  0.733333333333333 0.349019607843137])

yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% CHEN-DLT
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{CHEN_DLT_CONST}]
end
plot(yData,xData,'Color',[0.501960784313725  0.427450980392157 0.635294117647059],'LineWidth',2,'LineStyle','-','MarkerSize',6,'Marker','o','MarkerEdgeColor', [0.501960784313725  0.427450980392157 0.635294117647059])

yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% CHEN-FNS
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{CHEN_FNS_CONST}]
end
plot(yData,xData,'Color',[0.294117647058824  0.674509803921569 0.776470588235294],'LineWidth',2,'LineStyle','--','MarkerSize',6,'Marker','v','MarkerEdgeColor', [0.294117647058824  0.674509803921569 0.776470588235294])


yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% AML-DLT-2
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{AML_DLT_CONST}]
end
plot(yData,xData,'Color',[1 0 0],'LineWidth',2,'LineStyle','--','MarkerSize',6,'Marker','^','MarkerEdgeColor', [1 0 0])


yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% AML-FNS
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{AML_FNS_CONST}]
end
plot(yData,xData,'Color',[0 0.6 0.6],'LineWidth',2,'LineStyle','--','MarkerSize',6,'Marker','^','MarkerEdgeColor', [0 0.6 0.6])


yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% BA-DLT
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{BA_DLT_1_CONST}]
end
plot(yData,xData,'Color',[0 0.8 0],'LineWidth',2,'LineStyle','--','MarkerSize',6,'Marker','p','MarkerEdgeColor', [0 0.8 0])

yData = [];
xData = [];
for k = 1:numberOfNoiseLevels
% BA-FNS
yData = [yData k*noiseLevel];
xData = [xData listOfMeanErrorsForEachMethodAndSimulation{k}{BA_FNS_1_CONST}]
end
plot(yData,xData,'Color',[1 0.454901960784314 0],'LineWidth',2,'LineStyle','--','MarkerSize',6,'Marker','x','MarkerEdgeColor', [1 0.454901960784314 0])

legend('NDLT','FNS','CHEN-DLT','CHEN-FNS','AML-DLT','AML-FNS','BA-DLT','BA-FNS','Location','NorthWest');
xlabel('Noise Level');
ylabel('Mean RMS');
title('Symmetric Transfer Error From Ground Truth ');
saveas(myFig,'error_from_ground_truth_for_all_noise.fig')


end

