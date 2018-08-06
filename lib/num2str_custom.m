function [ str ] = num2str_custom( num )
%NUM2STR_CUSTOM Summary of this function goes here
%   Detailed explanation goes here
str = strrep(num2str(num),'.','_');

end

