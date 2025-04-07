function kernel = borderSupressionKernel(size, angle)
    % Create coordinate grid
    [X, Y] = meshgrid(linspace(-1, 1, size), linspace(1, -1, size));

    % Convert angle to radians
    theta = angle+pi/2;

    % Compute signed distance to the border line
    distance = X * cos(theta) + Y * sin(theta);

    % Normalize and create kernel with linear transition
    kernel = max(-1, min(1, distance));

    % Apply piecewise linear function
    kernel(kernel > 0) = 1 - kernel(kernel > 0);
    kernel(kernel < 0) = -1 - kernel(kernel < 0);

    gauss = gaussianFilter2D(size, size, 2, 2);
    kernel = gauss.*kernel;
end