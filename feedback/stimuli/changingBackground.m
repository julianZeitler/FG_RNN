function [array] = changingBackground(dim1, dim2, squareSize)
% Create an array with a linear transition from 0 to 1 from left to right
array = repmat(linspace(0, 1, dim1), dim2, 1);

% Calculate the starting and ending indices for the square
startIdx = (dim1 - squareSize) / 2 + 1;
endIdx = startIdx + squareSize - 1;

% Set the values inside the square to 0.5
array(startIdx:endIdx, startIdx:endIdx) = 0.5;
end

