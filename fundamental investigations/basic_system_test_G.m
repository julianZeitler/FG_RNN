path = 'output/feedback/20250611_Basic_System_Symmary/10';

iterations = 200;
alphaB = 0.36;
betaB = 1;
gammaB = 5;
zetaB = 0.9;
lambdaB = 5;
xiB = 0;

alphaG = 0.2;
betaG = 1;
gammaG = 2;
zetaG = 1;

tests = zeros(6, 2, 3);
%               GL Bias      GR Bias
tests(1,:,:) = [0   0   0  ; 0   0   0  ];
tests(2,:,:) = [0.5 0.5 0.5; 0.5 0.5 0.5];
tests(3,:,:) = [0   1   1  ; 0   0   0  ];
tests(4,:,:) = [0   1   1  ; 0   0.5 0.5];
tests(5,:,:) = [0   1   0  ; 0   0   1  ];
tests(6,:,:) = [0   0.5 0.5; 0   0   0.501  ];
tests = tests*0.1;

for test = 1:6
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
    
    figure;
    
    % Left subplot (full height)
    subplot(3,2,[1 3 5]); % spans all three rows in first column
    plot(1:iterations, BL_vals, 'b', 'DisplayName', 'BL');
    hold on;
    plot(1:iterations, BR_vals, 'r', 'DisplayName', 'BR');
    xlabel('Iteration');
    title('BL and BR over Iterations');
    legend('show', 'Location', 'southeast');
    grid on;
    
    % Right subplots (E, GL, GR)
    subplot(3,2,2);
    yyaxis left
    plot(1:iterations, E_vals, 'Color', [0 0.5 0], "DisplayName", "E");
    ylabel("E");
    ylim([0 1.2]);
    ax = gca;
    ax.YColor = [0 0.5 0];
    
    yyaxis right
    plot(1:iterations, GLBias_vals, 'b-', "DisplayName", "GLBias");
    hold on;
    plot(1:iterations, GRBias_vals, 'r-', "DisplayName", "GRBias");
    ylabel("GLBias/GRBias");
    ax = gca;
    ax.YColor = 'k';
    ylim([0, 1.2*max([max(GLBias_vals), max(GRBias_vals), 0.1])]);

    xlabel('Iteration');
    title('E, GLBias and GRBias over Iterations');
    legend('show', 'Location', 'southeast');
    grid on;
    
    subplot(3,2,[4, 6]);
    plot(1:iterations, GL_vals, 'b', 'DisplayName', 'GL');
    hold on;
    plot(1:iterations, GR_vals, 'r', 'DisplayName', 'GR');
    xlabel('Iteration');
    title('GL and GR over Iterations');
    legend('show', 'Location', 'southeast');
    grid on;


    % plot(1:iterations, GL_vals, 'g');
    % xlabel('Iteration');
    % ylabel('GL');lkjlkj;
    % title('GL over Iterations');
    % ylim([0 1.2]);
    % grid on;
    % 
    % subplot(3,2,6);
    % plot(1:iterations, GR_vals, 'b');
    % xlabel('Iteration');
    % ylabel('GR');
    % title('GR over Iterations');
    % ylim([0 1.2]);
    % grid on;
    
    set(gcf, "Position", [100, 100, 1200, 400]);
    
    % Add text box with parameter values
    paramText = sprintf('Parameters:\n\\alpha_B = %.2f\n\\beta_B = %.2f\n\\gamma_B = %.2f\n\\zeta_B = %.2f\n\\lambda_B = %.2f\n\\xi_B = %.2f\n\n\\alpha_G = %.2f\n\\beta_G = %.2f\n\\gamma_G = %.2f\n\\zeta_G = %.2f', ...
                        alphaB, betaB, gammaB, zetaB, lambdaB, xiB, alphaG, betaG, gammaG, zetaG);
    annotation(gcf, 'textbox', [0.01, 0.85, 0.2, 0.1], ...
               'String', paramText, ...
               'FitBoxToText', 'on', ...
               'BackgroundColor', 'white', ...
               'EdgeColor', 'black', ...
               'FontSize', 10);
    
    filename = sprintf('%d_BL_BR.png', test);
    saveas(gcf, fullfile(path, filename));
end

close all;