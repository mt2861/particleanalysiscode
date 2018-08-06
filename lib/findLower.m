function [lower_x, Npart_current] = findLower(lower_x_prev, NpartFunc, Ndesired, tol)
%FINDLOWER Summary of this function goes here
%   Detailed explanation goes here
lower_x = lower_x_prev;
Npart_current = NpartFunc(lower_x);            
diff_particles = Ndesired-Npart_current;
Nite = 0;
while abs(diff_particles)>tol
    lower_x = lower_x-sign(diff_particles);
    Npart_current = NpartFunc(lower_x);
    diff_particles = Ndesired-Npart_current;
    lower_x
    Nite = Nite+1;
    if Nite>1000
        disp('Ite to find lower_x')
        disp(Nite)
        tol = tol + 1;
    end
end
end

