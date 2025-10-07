function [params] = makeParams(options)
    arguments
        options.iterations = 20
        options.oris = deg2rad([0 22.5 45 67.5 90 112.5 135 157.5])
        options.num_scales = 3
        options.scale_step = 2
        options.R0 = 16
        options.filter_type {mustBeMember(options.filter_type, ["donut", "gauss"])} = "donut"
        options.G_alpha = 0.25
        options.G_beta = 1
        options.G_gamma = 0
        options.G_zeta = 0.9
        options.G_mu = 1
        options.G_nu = 0
        options.B_alpha = 0.36
        options.B_beta = 1
        options.B_gamma = 7
        options.B_zeta = 0.7
        options.B_lambda = 5
        options.B_mu = 0
        options.B_nu = 0
        options.B_integration_mode {mustBeMember(options.B_integration_mode, ["sum", "max"])} = "sum"
    end

    %% General parameters
    %number of iterations to run the model for
    params.iterations = options.iterations;
    params.oris = options.oris;
    params.num_ori = length(params.oris); % Oris are defined as the edge orientations
    params.num_scales = options.num_scales;
    params.scale_step = options.scale_step;
    params.R0 = options.R0;
    params.filter_type = options.filter_type;

    %% G-cell parameters
    params.G.alpha = options.G_alpha; % exp decay
    params.G.beta = options.G_beta; % scale
    params.G.gamma = options.G_gamma; % competition (not implemented)
    params.G.zeta = options.G_zeta; % shunting inhibition
    params.G.mu = options.G_mu; % inhibition from opposing B cells
    params.G.nu = options.G_nu; % scale competition

    %% B-Cell parameters
    params.B.alpha = options.B_alpha;
    params.B.beta = options.B_beta;
    params.B.gamma = options.B_gamma;
    params.B.zeta = options.B_zeta;
    params.B.lambda = options.B_lambda;
    params.B.mu = options.B_mu;
    params.B.nu = options.B_nu;
    params.B.integration_mode = options.B_integration_mode;

    %% Create filters
    for k = 0:params.num_scales-1
        GRF = makeGRF(params.R0*params.scale_step^k,params.oris+pi/2, params.filter_type);
        for ori = 1:length(params.oris)
            params.G.RF{k+1,ori} = GRF{ori};
            params.G.RF{k+1,ori+8} = GRF{ori+8};
        end
        params.G.inhibition_neighborhood{k+1} = gaussianFilter2D( ...
            3*params.R0*params.scale_step^(k+1)-mod(3*params.R0*params.scale_step^(k+1),2), ...
            3*params.R0*params.scale_step^(k+1)-mod(3*params.R0*params.scale_step^(k+1),2), ...
            1*params.R0*params.scale_step^(k+1), ...
            1*params.R0*params.scale_step^(k+1));
    end

    R1 = 1;
    params.B.FF.spatial_neighborhood_exc = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1, R1);
    params.B.FF.spatial_neighborhood_inh = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1*2, R1*2);
    params.B.FB.spatial_neighborhood = gaussianFilter2D(20*R1+(1-mod(20*R1,2)), 20*R1+(1-mod(20*R1,2)), R1*2, R1*2);
end

