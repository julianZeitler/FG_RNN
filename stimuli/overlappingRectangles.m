function [arr] = overlappingRectangles(dimX, dimY, height1, width1, height2, width2)
    arr = zeros(dimY, dimX);

    % First rectangle - centered, value 0.5
    center_y1 = round(dimY / 2);
    center_x1 = round(dimX / 2);
    y1_start = center_y1 - round(height1 / 2) + 1;
    y1_end = y1_start + height1 - 1;
    x1_start = center_x1 - round(width1 / 2) + 1;
    x1_end = x1_start + width1 - 1;
    arr(y1_start:y1_end, x1_start:x1_end) = 0.5;
    
    % Second rectangle - centered vertically, right-aligned with first, value 1
    center_y2 = round(dimY / 2);
    y2_start = center_y2 - round(height2 / 2) + 1;
    y2_end = y2_start + height2 - 1;
    x2_end = x1_end;  % Right edge aligns with first rectangle
    x2_start = x2_end - width2 + 1;
    arr(y2_start:y2_end, x2_start:x2_end) = 1;
end

