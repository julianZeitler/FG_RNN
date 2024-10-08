function [array] = verticalBarStimulus(dim1, dim2, width)
array = zeros(dim1, dim2);
array(:,floor(dim2/2-width/2):floor(dim2/2-width/2)+width) = 1;
end

