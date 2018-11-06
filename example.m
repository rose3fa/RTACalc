
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot the dispersion of two layer films, for the full matrix of different thicknesses for the two films
% date = '2018-11-05'; % Used for saving files
% suffix = 'Cu_WSe2_dispersionThicknessMatrix'; % Used to create new folders if repeating same runs on same day
% 
% lam0 = 1240/2.3;  % Smallest wavelength of interest (nm)
% lam1 = 1240/1.8;  % Largest wavelength of interest (nm)
% dlam = 10;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Ag' 'MoS2' 'air'};
% 
% film1thicknesses = [40 80];     % nm
% film2thicknesses = [1 5];   % nm
% 
% angle = 40:2:60;
% polarization = 1; 
% 
% single = [1 1];
% control = [1 0];
% 
% bilayerDispersionStudy(lam0, lam1, dlam, date, suffix, layers, film1thicknesses, film2thicknesses, angle, polarization, single, control);
% 
% clearvars date suffix lam0 lam1 dlam layers film1thicknesses film2thicknesses angle polarization control single;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot the dispersion of arbitrary layer films
% % Setup output directory
% date = '2018-11-05'; % Used for saving files
% suffix = 'Au_MoS2_dispersion'; % Used to create new folders if repeating same runs on same day
% experimentFolder = fullfile('Output', date, suffix);
% directory = fullfile(pwd, experimentFolder);
% mkdir(directory);
% 
% % Input parameters
% lam0 = 1240/2.5;  % Smallest wavelength of interest (nm)
% lam1 = 1240/1.5;  % Largest wavelength of interest (nm)
% dlam = 1;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Au' 'MoS2' 'air'};
% thicknesses = [50 5];     % nm
% angle = 1:1:89;
% polarization = 1; 
% dispersionX = 0;
% 
% % Used for AvsAngle
% wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
%                      % If lam0 + n*dlam doesn't equal lam1 for any n, it 
%                      % will round down to the nearest value to lam1
% 
% % Calculate the absorptivity vs angle
% [A, n_substrate]=AvsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl);
% 
% % Plot it
% plotDispersion(angle, wl, A, n_substrate, layers, thicknesses, polarization, directory, dispersionX);
% clearvars date suffix lam0 lam1 dlam layers angle polarization A directory dispersionX experimentFolder n_substrate note thicknesses wl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot R for various angles to help compare to spectrometer output
% % Setup output directory
% date = '2018-11-05'; % Used for saving files
% suffix = 'Au_MoS2'; % Used to create new folders if repeating same runs on same day
% experimentFolder = fullfile('Output', date, suffix);
% directory = fullfile(pwd, experimentFolder);
% mkdir(directory);
% 
% % Input parameters
% lam0 = 500;  % Smallest wavelength of interest (nm)
% lam1 = 850;  % Largest wavelength of interest (nm)
% dlam = 1;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Au' 'MoS2' 'air'};
% thicknesses = [35 5];     % nm
% angle = [42 47 52 60 70];
% polarization = 1; 
% 
% % Used for RvsAngle
% wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
%                      % If lam0 + n*dlam doesn't equal lam1 for any n, it 
%                      % will round down to the nearest value to lam1
% 
% % Calculate the reflectivity vs angle
% [R]=RvsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl);
% 
% % Plot it
% plotRvsAngle(wl, layers, lam0, lam1, angle, thicknesses, R, polarization, directory);
% clearvars date suffix lam0 lam1 dlam layers angle polarization R directory experimentFolder n_substrate note thicknesses wl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot A for various angles to help compare to spectrometer output
% % Setup output directory
% date = '2018-11-05'; % Used for saving files
% suffix = 'Au_MoS2'; % Used to create new folders if repeating same runs on same day
% experimentFolder = fullfile('Output', date, suffix);
% directory = fullfile(pwd, experimentFolder);
% mkdir(directory);
% 
% % Input parameters
% lam0 = 500;  % Smallest wavelength of interest (nm)
% lam1 = 850;  % Largest wavelength of interest (nm)
% dlam = 1;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Au' 'MoS2' 'air'};
% thicknesses = [35 5];     % nm
% angle = [42 47 52 60 65];
% polarization = 1; 
% 
% % Used for AvsAngle
% wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
%                      % If lam0 + n*dlam doesn't equal lam1 for any n, it 
%                      % will round down to the nearest value to lam1
% 
% % Calculate the A vs angle
% [A]=AvsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl);
% 
% % Plot it
% plotAvsAngle(wl, layers, lam0, lam1, angle, thicknesses, A, polarization, directory);
% clearvars date suffix lam0 lam1 dlam layers angle polarization A directory experimentFolder n_substrate note thicknesses wl;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot T for various angles to help compare to spectrometer output
% % Setup output directory
% date = '2018-11-05'; % Used for saving files
% suffix = 'Au_MoS2'; % Used to create new folders if repeating same runs on same day
% experimentFolder = fullfile('Output', date, suffix);
% directory = fullfile(pwd, experimentFolder);
% mkdir(directory);
% 
% % Input parameters
% lam0 = 500;  % Smallest wavelength of interest (nm)
% lam1 = 850;  % Largest wavelength of interest (nm)
% dlam = 1;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Au' 'MoS2' 'air'};
% thicknesses = [35 5];     % nm
% angle = [42 47 52 60 65];
% polarization = 1; 
% 
% % Used for TvsAngle
% wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
%                      % If lam0 + n*dlam doesn't equal lam1 for any n, it 
%                      % will round down to the nearest value to lam1
% 
% % Calculate the T vs angle
% [T]=TvsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl);
% 
% % Plot it
% plotTvsAngle(wl, layers, lam0, lam1, angle, thicknesses, T, polarization, directory);
% clearvars date suffix lam0 lam1 dlam layers angle polarization T directory experimentFolder n_substrate note thicknesses wl;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot Absorbance (-logT) for various angles to help compare to spectrometer output
% % Setup output directory
% date = '2018-11-05'; % Used for saving files
% suffix = 'Au_MoS2'; % Used to create new folders if repeating same runs on same day
% experimentFolder = fullfile('Output', date, suffix);
% directory = fullfile(pwd, experimentFolder);
% mkdir(directory);
% 
% % Input parameters
% lam0 = 500;  % Smallest wavelength of interest (nm)
% lam1 = 850;  % Largest wavelength of interest (nm)
% dlam = 1;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Au' 'MoS2' 'air'};
% thicknesses = [35 5];     % nm
% angle = [42 47 52 60 65];
% polarization = 1; 
% 
% % Used for AbsVsAngle
% wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
%                      % If lam0 + n*dlam doesn't equal lam1 for any n, it 
%                      % will round down to the nearest value to lam1
% 
% % Calculate the T vs angle
% [Abs]=AbsVsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl);
% 
% % Plot it
% plotAbsVsAngle(wl, layers, lam0, lam1, angle, thicknesses, Abs, polarization, directory);
% clearvars date suffix lam0 lam1 dlam layers angle polarization Abs directory experimentFolder n_substrate note thicknesses wl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Plot fitted dielectric functions
% % Use to verify interpolation/extrapolation is acceptable
% % Setup output directory
% date = '2018-11-05'; % Used for saving files
% suffix = 'Au_MoS2'; % Used to create new folders if repeating same runs on same day
% experimentFolder = fullfile('Output', date, suffix);
% directory = fullfile(pwd, experimentFolder);
% mkdir(directory);
% 
% % Input parameters
% lam0 = 500;  % Smallest wavelength of interest (nm)
% lam1 = 850;  % Largest wavelength of interest (nm)
% dlam = 1;         % Wavelength interval/resolution (nm)
% layers = {'fused silica' 'Au' 'MoS2' 'air'};
% thicknesses = [35 5];     % nm
% angle = 0;
% polarization = 1; 
% 
% % Used for plotFittedIndices
% wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
%                      % If lam0 + n*dlam doesn't equal lam1 for any n, it 
%                      % will round down to the nearest value to lam1
% 
% % Calculate the dielectric function fits
% [~, ~, ~, ~, ~, nData, nDataInterp, kData, kDataInterp, N] = thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle, polarization);
% % Plot it
% plotFittedIndices(wl, nData, nDataInterp, kData, kDataInterp, N, layers, lam0, lam1, directory);
% 
% clearvars date suffix lam0 lam1 dlam layers angle polarization Abs directory experimentFolder n_substrate note thicknesses wl kData kDataInterp N nData nDataInterp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot R, T, A
% Setup output directory
date = '2018-11-05'; % Used for saving files
suffix = 'Au_MoS2'; % Used to create new folders if repeating same runs on same day
experimentFolder = fullfile('Output', date, suffix);
directory = fullfile(pwd, experimentFolder);
mkdir(directory);

% Input parameters
lam0 = 500;  % Smallest wavelength of interest (nm)
lam1 = 850;  % Largest wavelength of interest (nm)
dlam = 1;         % Wavelength interval/resolution (nm)
layers = {'fused silica' 'Au' 'MoS2' 'air'};
thicknesses = [35 5];     % nm
angle = 60;
polarization = 1; 

% Used for AvsAngle
wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
                     % If lam0 + n*dlam doesn't equal lam1 for any n, it 
                     % will round down to the nearest value to lam1

% Calculate the RTA
[~, ~, R, T, A, ~, ~, ~, ~, ~] = thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle, polarization);

% Plot it
plotRTA(wl, layers, lam0, lam1, angle, thicknesses, R, T, A, polarization, directory)

clearvars date suffix lam0 lam1 dlam layers angle polarization R T A directory experimentFolder n_substrate note thicknesses wl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
