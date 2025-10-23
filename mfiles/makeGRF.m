function GRF = makeGRF(R0, oris, type)
%MAKEGRF Create RFs for G-cells
%   inputs:
%       R0: radius (mean of gaussian)
%       oris: orientations (only provide interval [0,pi), the rest is
%       inferred)
%       mode: "donut" or "gauss"
%
%   out:
%       RF: filter masks

[X, Y] = meshgrid(round(-2*R0):round(2*R0), round(2*R0):-1:round(-2*R0));
[phi, rho] = cart2pol(X, Y);

if type == "donut"
    values = normpdf(rho,R0,R0/3);
elseif type == "gauss"
    values = normpdf(rho,0,R0);
else
    error("Select either 'donut' or 'gauss' as mode.");
end

for ori=1:length(oris)
    sigma = pi/(length(oris));
    sigma2 = sigma^2;
    % Mask everything further than 2*sigma for performance improvement
    angle_mask_1 = abs(wrapToPi(phi-oris(ori)))<=2*sigma;
    angle_mask_2 = abs(wrapToPi(phi-oris(ori)+pi))-0.0001<=2*sigma;

    radial_gauss1 = 1/sqrt(2*pi*sigma2)*exp(-1/2*wrapToPi((phi-oris(ori))).^2/sigma2).*angle_mask_1;
    radial_gauss2 = 1/sqrt(2*pi*sigma2)*exp(-1/2*wrapToPi((phi-oris(ori)+pi)).^2/sigma2).*angle_mask_2;
    GRF{ori} = values.*radial_gauss1./sum(values.*radial_gauss1, "all");
    GRF{ori+length(oris)} = values.*radial_gauss2./sum(values.*radial_gauss2, "all");
end
end

