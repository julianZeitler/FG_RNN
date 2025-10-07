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
phi(ceil(size(phi, 1)/2), ceil(size(phi, 2)/2)) = NaN;

if type == "donut"
    values = normpdf(rho,R0,R0/3);
elseif type == "gauss"
    values = normpdf(rho,0,R0);
else
    error("Select either 'donut' or 'gauss' as mode.");
end

for ori=1:length(oris)
    angle_mask_1 = abs(wrapToPi(phi-oris(ori)))<=pi/(length(oris)*2);
    angle_mask_2 = abs(wrapToPi(phi-oris(ori)+pi))-0.0001<=pi/(length(oris)*2);
    GRF{ori} = values.*angle_mask_1./sum(values.*angle_mask_1, "all");
    GRF{ori+length(oris)} = values.*angle_mask_2./sum(values.*angle_mask_2, "all");
end
end

