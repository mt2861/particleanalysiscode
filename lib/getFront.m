function [ x_curr, indicies_ps_in_front ] = getFront(pos_COMs, num_ps_in_front, tol_front, step_move_front )
%GETFRONT Summary of this function goes here
%   Detailed explanation goes here
    indicies_ps_in_front = 1;
    x_curr = max(pos_COMs);
    % compares the actual number of particles in front with target
    % value num_particles_in_front
    diff_to_particle_front = (num_ps_in_front-numel(indicies_ps_in_front))/num_ps_in_front;
    while abs(diff_to_particle_front)>tol_front  
        x_curr = (1-sign(diff_to_particle_front)*step_move_front)*x_curr;
        %particles_front => indicies of particles in the top 50% of
        %the front
        indicies_ps_in_front = find(pos_COMs>x_curr);                   
        diff_to_particle_front = (num_ps_in_front-numel(indicies_ps_in_front))/num_ps_in_front;  
    end
    %xcurr now contains the xpos infront of which are 50% +/-
    %0.1% of the total number of particles

end

