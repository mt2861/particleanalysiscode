function [ centers ] = edges2centers( edges )
%EDGES Summary of this function goes here
%   Detailed explanation goes here
    centers = (edges(1:end-1)+edges(2:end))/2;
end

