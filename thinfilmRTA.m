function [wl, n_substrate, R, T, A, nData, nDataInterp, kData, kDataInterp, N] = thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle, polarization)


%% This script outputs spectral R, T, A for thin film stacks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%   About this script   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% How to use this function, and functions in general (basic usage)
    %      This file is a function, so you call it from either the command
    % window or another script with thinfilmRTA(...) and filling in each of
    % the variables according to the description below.
    %      To access the output variables either in the workspace or in
    % another script, you would use: [wl, n_substrate, R, T, A, nData,
    % nDataInterp, kData, kDataInterp, N] = thinfilmRTA(...), where each of
    % the variables in [] can be called anything (i.e., you don't have to
    % use the names above) and if a variable is unneeded it can be replaced
    % with ~ to save computation time. Further, to suppress output of each
    % variable in the Command Window, just add ; at the end of the command.

%% Input variables
    % lam0  Smallest wavelength of interest (nm)
    % lam1  Largest wavelength of interest (nm)
    % dlam  Wavelength interval/resolution (nm)
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
    % thicknesses is a vector (nm) that doesn't include the first or last
        % layers since they are always semi-infinite slabs
    % angle of incidence (deg)
    % polarization, 0 for TE (s-polarized), otherwise (any value other than
        % 0) TM (p-polarized)
    
%% Output variables
    % wl, a 1xn vector containing the wavelengths studied (nm)
    % n_substrate, a nx1 vector containing the real part of the index of
        % refraction of the substrate for each corresponding wavelength.
        % The underlying program neglects the imaginary part of the
        % substrate. Having n_substrate is using for calculating dispersion
        % plots
    % R, T, and A are the reflectance, trasmittance, and absorptivity (or
        % absorptance, but not absorbance...) at each wavelength, as a
        % decimal value, in a 1xn vector
    % nData and kData are cells (which are just arrays that can contain any
        % kind of information (not just numbers) with varying sizes) that
        % contain the actual imported real and imaginary refractive index
        % and wavelength data
    % nDataInterp and kDataInterp are the interpolated/extrapolated
        % refractive index data. Wavelengths are not included as they
        % correspond to the wl variable. Included so one can check that the
        % refractive index data being used makes sense, especially if it's
        % being extrapolated.
    % N is the number of layers in the stack
   
%% jreftran_rt.m
    % This script (thinfilmRTA.m) calls upon the jreftan_rt.m script
    % written by someone else
    %   function [r,t,R,T,A]=jreftran_rt(l,d,n,t0,polarization)
    %   l = free space wavelength, nm
    %   d = layer thickness vector, nm
    %       same length as n vector
    %       Ex: [NaN,200,NaN], 1st and last layers are semi-infinite slabs
    %   n = layer complex refractive index vector n = n + ik
    %       same length as d vector
    %       Ex: [1,0.9707 + 1.8562i,1], 1st and last layers are lossless
    %   t0 = angle of incidence, radians
    %   polarization should be 0 for TE (s-polarized), otherwise (any value
    	%   other than 0) TM (p-polarized)       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

%% Start importing and interpolating the dielectric functions at the
%  appropriate wavlengths
    % Note: dielectric function = refractive index

wl = lam0:dlam:lam1; % Creates a vector from the desired wavelength range.
                     % If lam0 + n*dlam doesn't equal lam1 for any n, it 
                     % will round down to the nearest value to lam1

% Infer file names from function input and import the Refractive Index data
delimiter = '\t';   % Dispersion files are tab delimited
N = length(layers); % Number of layers

% Preallocate memory for faster computation
    % Cells are used for generic arrays or vectors of strings (or any type
    % of data) while actual vectors or arrays can only contain numbers and
    % have to have a regular shape (e.g. cell {1 2 3; 4 5} is OK, but array
    % [1 2 3; 4 5] is not OK
    nFileName = cell(1,N);      % n for real part of refractive index
    kFileName = cell(1,N);      % k for imaginary part of refractive index
    nFullFileName = cell(1,N);
    kFullFileName = cell(1,N);
    nData = cell(1,N);
    kData = cell(1,N);
    nDataInterp = cell(1,N);
    kDataInterp = cell(1,N);


for k = 1:N
  nFileName{k} = strcat(layers{k},'_nm_n.txt');  % Create the n file names
  kFileName{k} = strcat(layers{k},'_nm_k.txt');  % Create the k file names
  % Create the n file names including the full path
  % Note that using fullfile to create directory paths should be OS
  % agnostic. E.g., '~/Refractive Indices/Ag_nm_n' would work on a Mac, but
  % not on a PC
  nFullFileName{k} = fullfile(pwd, 'Refractive Indices', nFileName{k});
  % Create the k file names including the full path
  kFullFileName{k} = fullfile(pwd, 'Refractive Indices', kFileName{k});
  
  % Import the data
      % Call wavelength (nm) as, eg,    nData{1}(:,1), and 
      % call index as, eg,              nData{1}(:,2)
      % where {1} says to use the first layer's data, : says to take all
      % rows, and 1 says to take the first column, so you'll end up with
      % nx1 vectors for this example
  nData{k} = importdata(nFullFileName{k}, delimiter);
  kData{k} = importdata(kFullFileName{k}, delimiter); 
  
  % Interpolate and extrapolate the data at the requested wavelengths
    % pchip will also be used to extrapolate. If extrapolation method isn't
    % good, add data points to dispersion text files to correct (if valid
    % to do so)
  nDataInterp{k} = interp1(nData{k}(:,1),nData{k}(:,2),wl','pchip');
  kDataInterp{k} = interp1(kData{k}(:,1),kData{k}(:,2),wl','pchip');
end

%% Construct a vector of n and k data 
%  Ignore wavelength data as the interpolated values all share the same
%  wavelength vector, wl
n = zeros(length(wl),N); % preallocate memory
n(:,1) = nDataInterp{1}; % First medium is lossless, so only taken n data
n(:,end) = nDataInterp{end}; % Final medium is lossless, so only taken n data
n_substrate = n(:,1); % Output substrate index to help with dispersion plots

% For other films, n = n + ik
for k = 2:N-1
n(:,k) = nDataInterp{k}+1i*kDataInterp{k};
end


%% Construct a vector of the thicknesses for each layer
d = zeros(1,length(thicknesses)+2); % preallocate memory
d(1) = NaN; % First and last media are semi-infinite
d(end) = NaN;

% Other films take thickness specified when calling function
for p = 2:length(d)-1
    d(p) = thicknesses(p-1);
end
 
%% Call the transfer matrix R, T, A calculator function and build spectra
t0 = pi*angle/180; % Convert angle to radians
l = length(wl);
% Preallocate memory
R = zeros(1,l);
T = zeros(1,l);
A = zeros(1,l);
for m = 1:l
       % Reflection and transmission coefficients, r, t, not used, so
       % replace output with ~. Can add back in if needed. 
       % [r(m),t(m),R(m),T(m),A(m)]=jreftran_rt(wl(m),d,n(m,:),t0,polarization);
           % full form of jreftran_rt
       [~,~,R(m),T(m),A(m)]=jreftran_rt(wl(m),d,n(m,:),t0,polarization);
end


