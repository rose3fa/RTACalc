function [Abs, n_substrate]=AbsVsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl)
% Preallocate memory
T = zeros(length(wl), length(angle));    
% Calculate A for every angle
for q = 1:length(angle) 
    [~, n_substrate,~,T(:,q),~, ~, ~, ~, ~, ~]=thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle(q), polarization);
end
Abs = -log10(T);