function filter = gaussianFilter2D(sizeX, sizeY, sigmaX, sigmaY)
    % Create a 2D Gaussian filter mask with parameterized size and standard
    % deviations. Values are scaled between 0 and 1
    %
    % sizeX - Width of the Gaussian filter (number of columns)
    % sizeY - Height of the Gaussian filter (number of rows)
    % sigmaX - Standard deviation in the X direction
    % sigmaY - Standard deviation in the Y direction
    
    % Create a grid of coordinates centered at zero
    [X, Y] = meshgrid(-(sizeX-1)/2:(sizeX-1)/2, -(sizeY-1)/2:(sizeY-1)/2);
    
    % Calculate the 2D Gaussian function
    filter = exp(- (X.^2 / (2 * sigmaX^2) + Y.^2 / (2 * sigmaY^2)));
    
    % Normalize the Gaussian filter so that its values sum to 1
    filter = filter/sum(filter, "all");
end