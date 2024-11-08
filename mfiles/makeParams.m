function [params] = makeParams()

%% General parameters
%number of iterations to run the model for
params.iterations = 5;

%ENTER EDGE MAP ORIENTATIONS (DEG) FOR 1ST QUADRANT ONLY INTO ORI
%   - computations are done on ori and ori + 90 degrees
ori = [0 22.5 45 67.5]; % 8 orientations
oris = deg2rad([ori (ori+90)]);

params.oris = oris;
params.num_ori = length(oris);

params.num_scales = 10;
params.scale_step = sqrt(2);

%% G-cell parameters
R0 = 2;
params.G.scale = 1;
params.G.exp_decay = 0.001;
params.G.inhibition_strength = 2;
params.G.inhibition_neighborhood{1} = gaussianFilter2D(10*R0-mod(10*R0,2), 10*R0-mod(10*R0,2), 3*R0, 3*R0);

% create RFs (Receptive Fields) for G-cells
dim = -3*R0:3*R0;
for ori = 1:params.num_ori
    [params.G.RF{1,ori}, params.G.RF{1,ori+8}] = makeVonMises(R0, params.oris(ori)+pi/2, dim, dim);
end
for k=2:params.num_scales
    for ori = 1:params.num_ori
        params.G.RF{k,ori} = imresize(params.G.RF{1, ori}, params.scale_step^(k-1), "cubic");
        params.G.RF{k,ori} = params.G.RF{k,ori}/sum(params.G.RF{k,ori}, "all");
        params.G.RF{k,ori+8} = imresize(params.G.RF{1, ori+8}, params.scale_step^(k-1), "cubic");
        params.G.RF{k,ori+8} = params.G.RF{k,ori+8}/sum(params.G.RF{k,ori+8}, "all");
    end

    params.G.inhibition_neighborhood{k} = gaussianFilter2D(10*R0*params.scale_step^k-mod(10*R0*params.scale_step^k,2), 10*R0*params.scale_step^k-mod(10*R0*params.scale_step^k,2), 3*R0*params.scale_step^k, 3*R0*params.scale_step^k);
end

%% B-cell parameters (FF=FeedForward, FB=FeedBack)
params.B.inhibition = 0.5;

R1 = 1;
params.B.exp_decay = 2;         % alpha
params.B.saturation = 0.5;      % zeta
params.B.FF.scale = 1;          % beta
params.B.FF.inhibition = 12;    % delta
params.B.FF.ori_norm = 1;
params.B.FF.spatial_neighborhood_exc = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1, R1);
params.B.FF.spatial_neighborhood_inh = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1*2, R1*2);

params.B.FB.scale = 2;          % lambda
params.B.FB.coarse_scale = 0.25;
params.B.FB.offset = 0.5;       % T_offset
params.B.FB.inhibition = 100;   % gamma
params.B.FB.spatial_neighborhood = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1*2, R1*2);

end

