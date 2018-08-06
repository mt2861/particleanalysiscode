function [ dat ] = readInputfile( filename_input )
%READINPUTFILE Summary of this function goes here
%   Detailed explanation goes here
fid  = fopen(filename_input,'r');
dat = struct;
while ~feof(fid)
    templine = fgetl(fid);
    while (isempty(templine) && ~feof(fid)) 
        templine = fgetl(fid);
    end
    if not(contains(templine,'#')) % if not commented
        line = strsplit(templine);
        if length(line) == 2 
            % 2 words => assign field
            if isnan(str2double(line{2}))
                dat.(line{1}) = line{2};
            else
                dat.(line{1}) = str2double(line{2});
            end
        else
            if contains(templine, 'periodic_length') 
                % of form "periodic_length [lx] [ly] [lz]"
                periodic_length = strsplit(templine);
                dat.lx = str2double(periodic_length{2});
                dat.ly = str2double(periodic_length{3});
                dat.lz = str2double(periodic_length{4});
            end
            if contains(templine, 'structure')
                % of form "structure [whoops{2}] [whoops{3}]"
                whoops = strsplit(templine);
                dat.structure = [whoops{2},' ',whoops{3}];
            end
        end
    end
end

end

