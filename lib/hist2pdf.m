function [ pdf, centers ] = hist2pdf( hist, edges )
%HIST2PDF Summary of this function goes here
%   Detailed explanation goes here
    centers = edges2centers(edges);
    pdf = hist/trapz(centers,hist);
end

