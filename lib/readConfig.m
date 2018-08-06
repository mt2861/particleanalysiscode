function [ data_traj ] = readConfig( file_blobs_tot, num_saved_steps, N )
%READCONFIG Summary of this function goes here
%   Detailed explanation goes here
    fid = fopen(file_blobs_tot,'r');
    delimiter = ' ';
    formatSpec = '%f%f%f%f%f%f%f%[\n\r]';
    dataArray = textscan( ...
        fid, ...
        formatSpec, ...
        'Delimiter',delimiter, ...
        'MultipleDelimsAsOne', true, ...
        'ReturnOnError', false, ...
        'EndOfLine', '\r\n' ...
        );
    data = [dataArray{1:end-1}];
    % every N+1th line is the number of particles so remove that line
    for j = 1:num_saved_steps
        data_traj(N*j-(N-1):N*j,:) = data((N+1)*j-(N-1):(N+1)*j,:);
    end 
end
