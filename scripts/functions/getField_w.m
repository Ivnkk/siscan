function [Spect,phase] = getField_w(scan)
%GETFIELD_W Summary of this function goes here
%   Detailed explanation goes here

Spect = abs(GGuess).^2./max(abs(GGuess).^2.);
phase = angle(scan);
phase = unwrap(phase);
end

