function [array] = vaseStimulus(dimX, dimY, min_width, max_width)
array = zeros(dimY, dimX);
amplitude = (max_width - min_width)/4;
offset = min_width/2 + amplitude;
for i=1:dimY
    low = round(dimX/2-(amplitude*sin(3*pi/dimY*i) + offset));
    high = round(dimX/2+(amplitude*sin(3*pi/dimY*i) + offset));
    array(i,low:high) = 1;
end
end

