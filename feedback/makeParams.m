function [params] = makeParams()
%MAKEPARAMS Summary of this function goes here
%   Detailed explanation goes here

%pyramid levels over which to compute
minLevel = 1;  %leave as 1
maxLevel = 10; %How many times do you want to downsample the pyramid
params.maxLevel = maxLevel;

%set downsampling mode
params.downSample = 'half'; %downsample by 'half' or 'full' octave. better performance at half

%number of iterations to run the model for
params.iterations = 10;

%intensity and color channel weights (80% intensity, 20% color)
params.gray = 0.8;
params.color = 0.2;

%ENTER EDGE MAP ORIENTATIONS (DEG) FOR 1ST QUADRANT ONLY INTO ORI
%   - computations are done on ori and ori + 90 degrees
ori = [0 22.5 45 67.5]; % 8 orientations
oris = deg2rad([ori (ori+90)]);

% G-cell parameters
params.G.inhibition = 0;
params.G.scale = 700;
params.G.minLevel = minLevel;
params.G.maxLevel = maxLevel;
params.G.oris = oris;
params.G.numOri = length(oris);
R0 = 25;

% create RF (Receptive Field) for G-cells
dim1 = -3*R0:3*R0;
dim2 = dim1;
for ori = 1:params.G.numOri
    [params.G.RF{ori}, params.G.RF{ori+8}] = makeVonMises(R0,params.G.oris(ori)+pi/2,dim1,dim2);
end

% B-cell parameters (FF=FeedForward, FB=FeedBack)
params.B.minLevel = minLevel;
params.B.maxLevel = maxLevel;
params.B.numOri = length(oris);
params.B.oris = oris;
params.B.inhibition = 0.5;

R1 = 1;
params.B.exp_decay = 2;      % alpha
params.B.saturation = 0.5;         % zeta
params.B.FF.scale = 1;          % beta
params.B.FF.inhibition = 12;     % delta
params.B.FF.spatial_neighborhood_exc = gaussianFilter2D(20*R1-1, 20*R1-1, R1, R1);
% params.B.FF.spatial_neighborhood_exc = gaussianFilter2D(3*R0, 3*R0, R0/8, R0/8);
% params.B.FF.spatial_neighborhood_exc = zeros(3, 3);
% params.B.FF.spatial_neighborhood_exc(2, 2) = 1;
params.B.FF.spatial_neighborhood_inh = gaussianFilter2D(20*R1-1, 20*R1-1, R1*2, R1*2);
% params.B.FF.spatial_neighborhood_inh = gaussianFilter2D(3*R0, 3*R0, R0/4, R0/4);

params.B.FB.scale = 2;           % lambda
params.B.FB.offset = 0.5;          % T_offset
params.B.FB.inhibition = 100;         % gamma
params.B.FB.spatial_neighborhood = gaussianFilter2D(20*R1-1, 20*R1-1, R1*2, R1*2);
% params.B.FB.spatial_neighborhood = gaussianFilter2D(3*R0, 3*R0, R0/4, R0/4);

end

