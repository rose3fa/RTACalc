function [A, n_substrate]=AvsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl)
% Preallocate memory
A = zeros(length(wl), length(angle));    
% Calculate A for every angle
for q = 1:length(angle) 
    [~, n_substrate,~,~,A(:,q), ~, ~, ~, ~, ~]=thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle(q), polarization);
end