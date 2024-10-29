function [RF1, RF2] = makeGRF(R0, theta0, sigma_phi)
%MAKEGRF Create RFs for G-cells using the Difference of Gaussians operator
% to create an off-center on-surround filter mask
%   inputs:
%       R0: radius of donut
%       theta0: direction of positive activation
%       phi: spanned angle
%       k: width of the donut
%
%   out:
%       RF1: filter mask
%       RF2: Filter mask offset by pi

[X, Y] = meshgrid(round(-2*R0):round(2*R0), round(2*R0):-1:round(-2*R0));
[phi, rho] = cart2pol(X, Y);

sigma_rho = R0^1.5/10;

phi = wrapTo2Pi(phi);

% Construct covariance matrix from eigenvectors and -values
Sigma = [sigma_phi 0; 0 sigma_rho];
mu1 = [theta0; R0];
mu2 = [wrapTo2Pi(theta0+pi); R0];

RF1 = zeros(size(phi));
RF2 = zeros(size(phi));

for i=1:size(phi,1)
    for j=1:size(phi,2)
        RF1(i,j) = 1/(2*pi)*1/sqrt(det(Sigma)) * exp(-1/2*([phi(i,j);rho(i,j)]-mu1).'*inv(Sigma)*([phi(i,j);rho(i,j)] - mu1));
        RF2(i,j) = 1/(2*pi)*1/sqrt(det(Sigma)) * exp(-1/2*([phi(i,j);rho(i,j)]-mu2).'*inv(Sigma)*([phi(i,j);rho(i,j)] - mu2));
    end
end

% Norm
RF1 = RF1/sum(RF1, "all");
RF2 = RF2/sum(RF2, "all");

end

