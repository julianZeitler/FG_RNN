function GRF = makeGRF(R0, oris)
%MAKEGRF Create RFs for G-cells
%   inputs:
%       R0: radius (mean of gaussian)
%       oris: orientations (only provide interval [0,pi), the rest is
%       inferred)
%
%   out:
%       RF: filter masks

[X, Y] = meshgrid(round(-2*R0):round(2*R0), round(2*R0):-1:round(-2*R0));
[phi, rho] = cart2pol(X, Y);

values = normpdf(rho,R0,R0/3);

for ori=1:length(oris)
    angle_mask_1 = abs(wrapToPi(phi-oris(ori)))<pi/16;
    angle_mask_2 = abs(wrapToPi(phi-oris(ori)+pi))<pi/16;
    GRF{ori} = values.*angle_mask_1./sum(values.*angle_mask_1, "all");
    GRF{ori+8} = values.*angle_mask_2./sum(values.*angle_mask_2, "all");
end
end

