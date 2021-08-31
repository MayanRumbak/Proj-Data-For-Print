function [R] = PNG2Yoav(PNG,SumNotation)
%this function convert PNG data into raster matrix notation by Yoav
% input: PNG as a PNG file with RGB values
%         SumNotation  is the convertion foramt
% output: R as a raster matrix
% we asume that thier is no 2 RGB with the same sum
PNG = sum(PNG,3);
R = zeros(size(PNG));
for i = 1:size(SumNotation(:,1))
    R = R + (PNG == SumNotation(i,1))*SumNotation(i,2);
end
end

%% example for code to creat the sum notation
% % the notation for the convertion from RGB to Yoav's
% Notation = [...
%     %R   G   B   notation
%     127 127 127 0 ;... % Air
%     240 240 240 1 ;... % White
%     0 90 158 2 ;... % Cyan
%     166 33 98 3 ;... % Magenta
%     200 189 3 4 ;... % Yellow
%     26 26 29 5 ;... % Black
%     0 0 0  6 ;... % Clear #########?????######
%     255 255 189 8 ;... % Support
%     ];
% % we are assuming that no 2 resins has the same sum...
% SumNotation(:,1) = sum(Notation(:,1:3),2);
% SumNotation(:,2) = Notation(:,4);