function [] = plotRTA(wl, layers, lam0, lam1, angle, thicknesses, R, T, A, polarization, directory)

%% Plot R, T, A
        % Create title text
        stack = [layers{1}, ' / '];
        note = [' ' layers{2} ': ' num2str(thicknesses(1)) ' nm'  ];
        for q=2:length(layers)-1
            stack = [stack,layers{q},' / '];
        end
        if length(layers) > 3
            for q = 3:length(layers)-1
                note = [note,', ' layers{q} ': ' num2str(thicknesses(q-1)) ' nm '  ];
            end
        end
        stack = [stack, layers{end}];
        % Polarization        
        if polarization==0
        pol =  'TE';
        else
        pol = 'TM';
        end
        % Angle of incidence
        aoi = num2str(angle);
        % Put it all together
        plotTitle = {['R, T, A at ' aoi ' deg with ' pol '-polarization'], ['Layers: ' stack], note, ' '};
       
        % Create name to save image
        saveStack = strcat(layers{1},'-');
        for q=2:length(layers)-1
            saveStack = strcat(saveStack,layers{q},'-');
        end
        saveStack = strcat(saveStack, layers{end});
        saveTitle = ['R, T, A at ' aoi ' deg with ' pol '-polarization ', saveStack, ', ' note];
        saveTitle = replace(saveTitle, ':', '');  % colons aren't allowed in Windows filenames, so remove them

    % Set up plot
    font = 24;
    xLabel = 'Wavelength (nm)';
    yLabel = 'R, T, or A';
    dataLabel = {'R' 'T' 'A'};

    Plot = figure;
    set(Plot, 'Position', [1 1 1400 860]);
    axes('FontSize', font)
    xlabel(xLabel, 'FontSize', font)
    ylabel(yLabel, 'FontSize', font)
    
    % Plot data
    hold on
    a = scatter(wl, R, 20, [0.5 0.5 0.5]);
    b = scatter(wl, T, 20, [0.722 0.451 0.2]);
    c = scatter(wl, A, 100, 'k', 'filled');
 
    legend(dataLabel,'location', 'northeast')
    title(plotTitle, 'FontSize', font+2)
    axis([lam0 lam1 -inf inf])
    hold off
    
    saveas(Plot, fullfile(directory,[saveTitle '.png']));
