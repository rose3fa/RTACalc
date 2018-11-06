function [R]=RvsAngle(lam0, lam1, dlam, layers, thicknesses, angle, polarization, wl)
% Preallocate memory
R = zeros(length(wl), length(angle));    
% Calculate A for every angle
for q = 1:length(angle) 
    [~,~,R(:,q),~,~, ~, ~, ~, ~, ~]=thinfilmRTA(lam0, lam1, dlam, layers, thicknesses, angle(q), polarization);
end