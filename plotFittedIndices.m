function [] = plotFittedIndices(wl, nData, nDataInterp, kData, kDataInterp, N, layers, lam0, lam1, directory)

%%  Use to check validity of interpolation/extrapolation
% Setup
font = 24;
xLabel = 'Wavelength (nm)';
yLabelL = 'n';
yLabelR = 'k';
dataLabel = {'n interpolated' 'n data' 'k interpolated' 'k data' };
for k = 1:N
Plot = figure;
set(Plot, 'Position', [1 1 1400 860]);
axes('FontSize', font)
xlabel(xLabel, 'FontSize', font)
ylabel(yLabelL, 'FontSize', font)
hold on
% Plot n data
scatter(wl', nDataInterp{k}, 200, [1 165/255 0])
scatter(nData{k}(:,1), nData{k}(:,2), 200,'b', 'filled')
% Plot k data
yyaxis right
ylabel(yLabelR, 'FontSize', font)
scatter(wl', kDataInterp{k}, 200, 'r')
scatter(kData{k}(:,1), kData{k}(:,2), 200, 'g', 'filled')
legend(dataLabel,'location', 'northeast')
title(layers{k}, 'FontSize', font+2)
axis([lam0 lam1 -inf inf])
hold off
plotTitle = strcat('Interpolated dielectric function of_', layers{k});
saveas(Plot, fullfile(directory,[plotTitle '.png']));
end