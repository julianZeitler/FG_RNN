function [array] = vaseStimulus(dim1, dim2, min_width, max_width)
array = zeros(dim1, dim2);
amplitude = (max_width - min_width)/4;
offset = min_width/2 + amplitude;
for i=1:dim1
    low = round(dim2/2-(amplitude*sin(3*pi/dim1*i) + offset));
    high = round(dim2/2+(amplitude*sin(3*pi/dim1*i) + offset));
    array(i,low:high) = 1;
end
end

