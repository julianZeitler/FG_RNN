function [params] = makeParams()

%% General parameters
%number of iterations to run the model for
params.iterations = 5;

%ENTER EDGE MAP ORIENTATIONS (DEG) FOR 1ST QUADRANT ONLY INTO ORI
%   - computations are done on ori and ori + 90 degrees
ori = [0 22.5 45 67.5]; % 8 orientations
oris = deg2rad([ori (ori+90)]);

params.oris = oris;
params.num_ori = length(oris); % Oris are defined as the edge orientations

params.num_scales = 6;
params.scale_step = 2;
params.R0 = 2;

%% G-cell parameters
params.G.scale = 0.3;
params.G.exp_decay = 0.001;
params.G.inhibition_strength = 2;
params.G.inhibitory_input_strength = 0.2;

% create RFs (Receptive Fields) for G-cells
for k = 0:params.num_scales-1
    GRF = makeGRF(params.R0*params.scale_step^k,oris+pi/2);
    for ori = 1:length(oris)
        params.G.RF{k+1,ori} = GRF{ori};
        params.G.RF{k+1,ori+8} = GRF{ori+8};
    end
    params.G.inhibition_neighborhood{k+1} = gaussianFilter2D( ...
        10*params.R0*params.scale_step^(k+1)-mod(10*params.R0*params.scale_step^(k+1),2), ...
        10*params.R0*params.scale_step^(k+1)-mod(10*params.R0*params.scale_step^(k+1),2), ...
        3*params.R0*params.scale_step^(k+1), ...
        3*params.R0*params.scale_step^(k+1));
end

%% B-cell parameters (FF=FeedForward, FB=FeedBack)
params.B.inhibition = 0.5;

R1 = 1;
params.B.exp_decay = 2;         % alpha
params.B.saturation = 0.5;      % zeta
params.B.FF.scale = 1;          % beta
params.B.FF.inhibition = 7;    % delta
params.B.FF.ori_norm = 1;
params.B.FF.spatial_neighborhood_exc = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1, R1);
params.B.FF.spatial_neighborhood_inh = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1*2, R1*2);

params.B.FB.scale = 2;          % lambda
params.B.FB.coarse_scale = 0.25;
params.B.FB.offset = 0.5;       % T_offset
params.B.FB.inhibition = 300;   % gamma
params.B.FB.spatial_neighborhood = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1*2, R1*2);
params.B.FB.border_supression = 15;

end

