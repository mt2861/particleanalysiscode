function [ num ] = str2num_custom( str )
%STR2NUM_CUSTOM Summary of this function goes here
%   Detailed explanation goes here
num = str2num(strrep(str,'_','.'));

end

