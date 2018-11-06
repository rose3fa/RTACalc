function [Ac, Aavg, n_substrate]=aveOfFilms(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl)
%% This script calculates the dispersion of an arbitrary thin film stack and the dispersion of the average of each individual film
% Stacks can have any number of layers
    %% Inputs
    % lam0 Smallest wavelength of interest (nm)
    % lam1 Largest wavelength of interest (nm)
    % dlam Wavelength interval/resolution (nm)
    %                      
    % layers, e.g. {'fused silica' 'Ag' 'MoS2' 'air'};   
        % call layers{1} to get the string  
        % Must match first part of names in Dispersions folder
    % thicknesses (nm), e.g. [55 0.65]
        % Not including the substrate and top layer (which are assumed to
        % be semi-infinite)
    % angle (deg), e.g. angle = 0.1:0.1:89.9
        % A vector
        % This is the set of angles used to create the dispersion plot
        % These angles are the angles of incidence inside the substrate
    % polarization 
        % 1 for p-polarized
        % 0 for s-polarized
    % wl is just lam0:dlam:lam1, may be faster to pass it since it's
    % already calculated 
        
    %% Outputs
    % Ac is matrix of full stack absorptivity with dimensions wl x angle
    % Aavg is matrix of control stack absorptivity with dimensions wl x angle
    % n_substrate
        
        
%% How the calculation is done
% Create a 3 layer structure where the first and last layers are the same
% as the stack under study. The middle layer will be one of the middle
% layers in the stack under study. This 3-layer stack will be a "control"
% stack. All possible control stacks will be calculated, i.e. for all middle
% films of device under study. Then the control stacks will be averaged to
% serve as a "control" against the device under study. E.g., 2D TMDs on
% plasmonic metal exhibit strong coupling, so the "real" device will show
% much higher absorption in regions then the control device. 


%% Set up calculation  
% The first and last layers of the controls are the same as the "real"
% stack
controlLayers(1) = layers(1);
controlLayers(3) = layers(end);

% Preallocate memory
A = zeros(length(wl), length(angle), length(layers));  % Controls
Ac = zeros(length(wl), length(angle));                 % Full stack

for q = 1:length(angle) 
    [~, ~,~,~,Ac(:,q), ~, ~, ~, ~, ~]=thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle(q), polarization);
    for m = 1:(length(layers)-2)
        controlLayers(2) = layers(m+1); 
        controlThickness = thicknesses(m);
        [~,n_substrate,~,~,A(:,q,m),~,~,~,~,~]=thinfilmRTA(lam0, lam1, dlam, controlLayers, controlThickness, angle(q), polarization);
    end
end
Aavg= sum(A,3)./(length(layers)-2);


