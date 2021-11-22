function [Length] = Mm2DPI(Length,DPI)
%convert length in mm into picxels using DPI and inchs
% Mm2DPI(Length,DPI)
Length = round(Length/25.4*DPI);
end

