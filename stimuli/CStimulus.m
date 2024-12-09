function [arr] = CStimulus(dim, CSize, CWidth, amplitude)
%CSTIMULUS Summary of this function goes here
%   Detailed explanation goes here

arr = zeros(dim, dim);

startIdx = round((dim - CSize) / 2);
endIdx = round(startIdx + CSize);

arr(startIdx:endIdx, startIdx:endIdx) = amplitude;
arr(startIdx+CWidth:endIdx-CWidth,startIdx+CWidth:endIdx) = 0;

end

