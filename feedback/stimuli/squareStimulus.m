function [arr] = squareStimulus(dim1,dim2,squareSize)
%SQUARE Summary of this function goes here
%   Detailed explanation goes here
arr = zeros(dim1, dim2);

% Calculate the starting and ending indices for the square
startIdx = (dim1 - squareSize) / 2 + 1;
endIdx = startIdx + squareSize - 1;

arr(startIdx:endIdx, startIdx:endIdx) = 1;
end

