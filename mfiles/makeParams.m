function [params] = makeParams()

%% General parameters
%number of iterations to run the model for
params.iterations = 10;

%ENTER EDGE MAP ORIENTATIONS (DEG) FOR 1ST QUADRANT ONLY INTO ORI
%   - computations are done on ori and ori + 90 degrees
ori = [0 22.5 45 67.5]; % 8 orientations
oris = deg2rad([ori (ori+90)]);

params.oris = oris;
params.num_ori = length(oris);

%% G-cell parameters
R0 = 2;
params.G.num_scales = 10;
params.G.scale_step = sqrt(2);
params.G.scale = 1;
params.G.exp_decay = 0.0015;
params.G.inhibition_strength = 0.3;
params.G.inhibition_neighborhood = gaussianFilter2D(10*R0-mod(10*R0,2), 10*R0-mod(10*R0,2), 3*R0, 3*R0);

% create RFs (Receptive Fields) for G-cells
for k=1:params.G.num_scales
    R = R0*params.G.scale_step^(k-1);
    dim1 = round(-3*R:3*R);
    dim2 = dim1;
    for ori = 1:params.num_ori
        [params.G.RF{k,ori}, params.G.RF{k,ori+8}] = makeVonMises(R,params.oris(ori)+pi/2,dim1,dim2);
    end
end

%% B-cell parameters (FF=FeedForward, FB=FeedBack)
params.B.inhibition = 0.5;

R1 = 1;
params.B.exp_decay = 2;         % alpha
params.B.saturation = 0.5;      % zeta
params.B.FF.scale = 1;          % beta
params.B.FF.inhibition = 12;    % delta
params.B.FF.ori_norm = 1;
params.B.FF.spatial_neighborhood_exc = gaussianFilter2D(20*R1-mod(20*R1,2), 20*R1-mod(20*R1,2), R1, R1);
params.B.FF.spatial_neighborhood_inh = gaussianFilter2D(20*R1-mod(20*R1,2), 20*R1-mod(20*R1,2), R1*2, R1*2);

params.B.FB.scale = 2;          % lambda
params.B.FB.offset = 0.5;       % T_offset
params.B.FB.inhibition = 100;   % gamma
params.B.FB.spatial_neighborhood = gaussianFilter2D(20*R1-mod(20*R1,2), 20*R1-mod(20*R1,2), R1*2, R1*2);

end

