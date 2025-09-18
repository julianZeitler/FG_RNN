path = 'output/feedback/20250611_Basic_System_Symmary/12';

iterations = 200;
alphaB = 0.36;
betaB = 1;
gammaB = 5;
zetaB = 1.4;
lambdaB = 5;
xiB = 0;

alphaG = 0.2;
betaG = 1;
gammaG = 2;
zetaG = 1;

count = 1;
for l_bias = 0:0.05:1
    for r_bias = 0:0.05:1
        tests(count,:,:) = [0 l_bias l_bias; 0 0 l_bias+r_bias];
        count = count + 1;
    end
end

tests = tests*0.1;

delta_v1 = zeros(1, count-1);
delta_bias2 = zeros(1, count-1);
switched = zeros(1, count-1);

for test = 1:(count-1)
    E = 0;
    GL = 0;
    GR = 0;
    BL = 0;
    BR = 0;

    GLBias = 0;
    GRBias = 0;
    
    BL_vals = zeros(1, iterations);
    BR_vals = zeros(1, iterations);
    E_vals = zeros(1, iterations);
    GL_vals = zeros(1, iterations);
    GR_vals = zeros(1, iterations);
    GLBias_vals = zeros(1, iterations);
    GRBias_vals = zeros(1, iterations);

    for i = 1:iterations
        if i == 2
            E = 1;
            GLBias = tests(test,1,1);
            GRBias = tests(test,2,1);
        elseif i == 10
            GLBias = tests(test,1,2);
            GRBias = tests(test,2,2);
        elseif i == 50
            GLBias = tests(test,1,3);
            GRBias = tests(test,2,3);
        end
        BL_prev = BL;
        BR_prev = BR;
        GL_prev = GL;
        GR_prev = GR;
        % gammaB*BR_prev*(1 + xiB * GR_prev) + ...
        % gammaB*(BR_prev + xiB * GR_prev) + ...
        BL = betaB*E*(1+lambdaB*GL_prev)/( ...
            alphaB + ...
            gammaB*(BR_prev + xiB * GR_prev) + ...
            zetaB*E*(1+lambdaB*GL_prev) ...
        );
        % gammaB*BL_prev*(1 + xiB * GL_prev) + ...
        % gammaB*(BL_prev + xiB * GL_prev) + ...
        BR = betaB*E*(1+lambdaB*GR_prev)/( ...
            alphaB + ...
            gammaB*(BL_prev + xiB * GL_prev) + ...
            zetaB*E*(1+lambdaB*GR_prev) ...
        );

        GL = betaG*(BL+GLBias)/(alphaG + gammaG*GR_prev + zetaG*(BL+GLBias));
        GR = betaG*(BR+GRBias)/(alphaG + gammaG*GL_prev + zetaG*(BR+GRBias));
    
        BL_vals(i) = BL;
        BR_vals(i) = BR;
        E_vals(i) = E;
        GL_vals(i) = GL;
        GR_vals(i) = GR;
        GLBias_vals(i) = GLBias;
        GRBias_vals(i) = GRBias;
    end
    delta_v1(test) = abs(BL_vals(49) - BR_vals(49));
    delta_bias2(test) = tests(test,2,3) - tests(test,1,3);
    if (GR_vals(iterations) > GL_vals(iterations))
        switched(test) = 1;
    else
        switched(test) = 0;
    end
end

colors = [0 0 1;   % blue for not switched (0)
          1 0 0];  % red for switched (1)

% Map each point's color using switched values
point_colors = colors(switched + 1, :);

figure;
scatter(delta_v1, delta_bias2, 40, point_colors, 'filled');

xlabel('\Delta V_1');
ylabel('\Delta bias_2');
title('Switching Behavior');
grid on;

saveas(gcf, fullfile(path, "switch_space.png"));
close all;