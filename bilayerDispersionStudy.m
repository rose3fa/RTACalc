function [] =  bilayerDispersionStudy(lam0, lam1, dlam, date, suffix, layers, film1thicknesses, film2thicknesses, angle, polarization, single, control)
%% Plots the dispersion of thin film stacks with 2 films sandwiched between 2 semi-infinite media.

%% Input variables
    % lam0  Smallest wavelength of interest (nm)
    % lam1  Largest wavelength of interest (nm)
    % dlam  Wavelength interval/resolution (nm)
    % date
    % suffix    Outputs will be saved in Output/date/suffix
    % layers is a cell array entered e.g. {'fused silica' 'Ag' 'air'}
        % The first and last layers will be modeled as lossless
            % semi-infinite slabs; they can have dispersion but only the
            % real part of the refractive index will be taken.
        % Light is incident from first layer
        % Dielectric function are pulled from the 'Refractive Indices'
            % folder in the same directory as this function.
        % Naming of refractive index files is, e.g. 'Ag_nm_n' where all
            % files use 'nm,' 'n' can be 'n' or 'k,' and the material name
            % changes.
        % Call the layer based on the name of the file before the first
            % underscore, e.g. 'Ag'
    % film1thicknesses is a vector (nm) that includes the thicknesses you
        % wish to study for layer 2
    % film2thicknesses is a vector (nm) that includes the thicknesses you
        % wish to study for layer 3
    % angle     angles of incidence (deg), a vector to set the bounds of
        % dispersion plot
    % polarization, 0 for TE (s-polarized), otherwise (any value other than
        % 0) TM (p-polarized)
    % single     form: [(plot it?) (theta or kx?)]
               % plot it? 1 to plot the dispersion of the stack, 0 to skip
               % kx or deg? 1 for theta or 0 for kx for dispersion x axis    
    % control    form: [(plot it?) (theta or kx?)]
               % plot it? 1 to plot the dispersion of the stack vs the control, 0 to skip
               % kx or deg? 1 for theta or 0 for kx for dispersion x axis

% Makes calls to aveOfFilms and plot2dispersions
    
% Setup 
experimentFolder = fullfile('Output', date, suffix);
wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
                     % If lam0 + n*dlam doesn't equal lam1 for any n, it 
                     % will round down to the nearest value to lam1
directory = fullfile(pwd, experimentFolder);
mkdir(directory);

% Plot for all combinations of film thicknesses
for m = 1:length(film1thicknesses)
    thicknesses(1) =film1thicknesses(m);
    for k = 1:length(film2thicknesses)
        thicknesses(2) = film2thicknesses(k);
        % Get the stack and control absorptivity
        [Ac, Aavg, n_substrate] = aveOfFilms(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl);
        %% Plot it
        if single(1)==1
        plotDispersion(angle, wl, Ac, n_substrate, layers, thicknesses, polarization, directory, single(2));
        else
        end
        if control(1)==1
        plot2dispersions(angle, wl, Ac, Aavg, n_substrate, layers, thicknesses, polarization, directory, control(2));
        else
        end
    end
end