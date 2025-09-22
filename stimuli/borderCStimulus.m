function [arr] = borderCStimulus(dimX, dimY, CSizeX, CSizeY, CWidth, amplitude)
    arr = zeros(dimY, dimX);
    corner0 = [round((dimY - CSizeY) / 2), round((dimX - CSizeX) / 2)];
    corner1 = [corner0(1), corner0(2) + CSizeX];
    corner2 = [corner1(1)+CWidth, corner1(2)];
    corner3 = [corner2(1), corner2(2)-CSizeX+CWidth];
    corner4 = [corner3(1)+CSizeY-2*CWidth, corner3(2)];
    corner5 = [corner4(1), corner4(2)+CSizeX-CWidth];
    corner6 = [corner5(1)+CWidth, corner5(2)];
    corner7 = [corner6(1), corner6(2)-CSizeX];
    

    arr(corner1(1), corner0(2):corner1(2)) = amplitude;
    arr(corner1(1):corner2(1), corner2(2)) = amplitude;
    arr(corner3(1), corner3(2):corner2(2)) = amplitude;
    arr(corner3(1):corner4(1), corner4(2)) = amplitude;
    arr(corner5(1), corner4(2):corner5(2)) = amplitude;
    arr(corner5(1):corner6(1), corner6(2)) = amplitude;
    arr(corner7(1), corner7(2):corner6(2)) = amplitude;
    arr(corner0(1):corner7(1), corner0(2)) = amplitude;
end
