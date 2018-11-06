function [] = plotAbsVsAngle(wl, layers, lam0, lam1, angle,thicknesses, Abs, polarization, directory)

%% Plot Abs(theta)[OD]
% Absorbance is -log10(T), not the same as absorptivity (aka absorptance)
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
        % Put it all together
        plotTitle = {['Absorbance with ' pol '-polarization'], ['Layers: ' stack], note, ' '};
       
        % Create name to save image
        saveStack = strcat(layers{1},'-');
        for q=2:length(layers)-1
            saveStack = strcat(saveStack,layers{q},'-');
        end
        saveStack = strcat(saveStack, layers{end});
        saveTitle = ['Absorbance with ' pol '-polarization ', saveStack, ', ' note];
        saveTitle = replace(saveTitle, ':', '');  % colons aren't allowed in Windows filenames, so remove them   

    % Set up plot
    font = 24;
    xLabel = 'Wavelength (nm)';
    yLabel = 'Absorbance (OD)';
    for l = 1:length(angle)
        aoi{l} = num2str(angle(l));
    end
    
    Plot = figure;
    set(Plot, 'Position', [1 1 1400 860]);
    axes('FontSize', font)
    xlabel(xLabel, 'FontSize', font)
    ylabel(yLabel, 'FontSize', font)
    
    % Plot data
    hold on
    for q = 1:length(angle)
        scatter(wl, Abs(:,q),(100+20*(q-1)), 'filled');
    end
 
    a=legend(aoi,'location', 'northeast');
    a.Title.String = 'Angle of Incidence (deg)';
    title(plotTitle, 'FontSize', font+2)
    axis([lam0 lam1 -inf inf])
    hold off
    
    saveas(Plot, fullfile(directory,[saveTitle '.png']));

