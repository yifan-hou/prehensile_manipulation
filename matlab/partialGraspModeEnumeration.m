function [e_modes, h_modes] = partialGraspModeEnumeration(CP_W_e, CN_W_e, CP_H_h, CN_H_h)

% Enumerate contact modes
e_modes = contact_mode_enumeration(CP_W_e(1:2,:), CN_W_e(1:2,:), true);
h_modes = contact_mode_enumeration(CP_H_h(1:2,:), CN_H_h(1:2,:), true);

% get rid of all separation modes
e_modes(:, all(e_modes == 0, 1)) = [];
h_modes(:, all(h_modes == 0, 1)) = [];

% add all fixed mode if necessary
if sum(all(e_modes == 1, 1)) == 0
    e_modes = [e_modes ones(size(e_modes,1),1)];
end
if sum(all(h_modes == 1, 1)) == 0
    h_modes = [h_modes ones(size(h_modes,1),1)];
end

disp(['Enumerated modes: ' num2str(size(e_modes, 2)*size(h_modes, 2)) ' modes.']);
