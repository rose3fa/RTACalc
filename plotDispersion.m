function [] = plotDispersion(angle, wl, A, n_substrate, layers, thicknesses, polarization, directory, dispersionX)
%% This function plots dispersion of a stack
    % dispersionX, 0 for kx or 1 for theta for plotting dispersion as energy vs kx
        % or angle

%% Compute E and kx based on Abs, wl, and theta
Abs = A;
eV = 1240 ./ wl';   % Convert wavelengths to energy (eV)
kx = 2*pi*sin(pi*angle/180).*n_substrate./wl';
    
%% Setup titles
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
        plotTitle = {['Dispersion with ' pol '-polarization'], ['Layers: ' stack], note, ' '};
       
        % Create name to save image
        saveStack = strcat(layers{1},'-');
        for q=2:length(layers)-1
            saveStack = strcat(saveStack,layers{q},'-');
        end
        saveStack = strcat(saveStack, layers{end});
        saveTitle = ['Dispersion with ' pol '-polarization ', saveStack, ', ' note];
        saveTitle = replace(saveTitle, ':', '');  % colons aren't allowed in Windows filenames, so remove them
    
    %% Set up plot
    set(0,'DefaultFigureVisible','off'); % Don't display the plot--just save it
    font = 24;
    
    % Plot vs kx or angle
    if dispersionX == 0
    xLabel = '\textbf{Wavenumber} $\frac{2\pi}{\lambda}n\sin{\theta} \ \  (nm^{-1})$';
    xAxis = kx;
    else 
    xLabel = '\textbf{Angle} $\theta \ \  (degrees)$'; 
    xAxis = angle;
    end
    yLabel = '\textbf{Photon Energy} $(eV)$';

   
    Plot = figure;
    set(Plot, 'Position', [1 1 1400 860]);
    axes('FontSize', font)  ;
    xlabel(xLabel, 'FontSize', font, 'Interpreter','latex');
    ylabel(yLabel, 'FontSize', font, 'Interpreter','latex');
    
hold on
    plot = pcolor(xAxis, eV, Abs);
    plot.EdgeColor = 'none';
    colormap('hot')
    c=colorbar;
    c.Label.String = 'Absorptivity';
    title(plotTitle, 'FontSize', font+1)
hold off
    saveas(Plot, fullfile(directory,[saveTitle '.png']));
