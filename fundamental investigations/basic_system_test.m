path = 'output/feedback/20250517_Shunting_Inhibition';

iterations = 30;
alpha = 1;
beta = 1;
gamma = 1;
zeta = 1;

tests = zeros(6, 2, 3);
%               FBL          FBR
tests(1,:,:) = [0   0   0  ; 0   0   0  ];
tests(2,:,:) = [0.5 0.5 0.5; 0.5 0.5 0.5];
tests(3,:,:) = [0   1   1  ; 0   0   0  ];
tests(4,:,:) = [0   1   1  ; 0   0.5 0.5];
tests(5,:,:) = [0   1   0  ; 0   0   1  ];
tests(6,:,:) = [0   0.5 0.5; 0   0   1  ];

for test = 1:6
    E = 0;
    FBL = 0;
    FBR = 0;
    BL = 0;
    BR = 0;
    
    BL_vals = zeros(1, iterations);
    BR_vals = zeros(1, iterations);
    E_vals = zeros(1, iterations);
    FBL_vals = zeros(1, iterations);
    FBR_vals = zeros(1, iterations);

    for i = 1:iterations
        if i == 2
            E = 1;
            FBL = tests(test,1,1);
            FBR = tests(test,2,1);
        elseif i == 10
            FBL = tests(test,1,2);
            FBR = tests(test,2,2);
        elseif i == 15
            FBL = tests(test,1,3);
            FBR = tests(test,2,3);
        end
        BL_prev = BL;
        BR_prev = BR;
    
        BL = beta*E*(1+FBL)/(alpha + zeta*E*(1+FBL) + gamma*BR_prev);
        BR = beta*E*(1+FBR)/(alpha + zeta*E*(1+FBR) + gamma*BL_prev);
    
        BL_vals(i) = BL;
        BR_vals(i) = BR;
        E_vals(i) = E;
        FBL_vals(i) = FBL;
        FBR_vals(i) = FBR;
    end
    
    figure;
    
    % Left subplot (full height)
    subplot(3,2,[1 3 5]); % spans all three rows in first column
    plot(1:iterations, BL_vals, 'DisplayName', 'BL');
    hold on;
    plot(1:iterations, BR_vals, 'DisplayName', 'BR');
    xlabel('Iteration');
    ylabel('Value');
    title('BL and BR over Iterations');
    legend('show', 'Location', 'southeast');
    grid on;
    
    % Right subplots (E, FBL, FBR)
    subplot(3,2,2);
    plot(1:iterations, E_vals, 'r');
    xlabel('Iteration');
    ylabel('E');
    title('E over Iterations');
    ylim([0 1.2]);
    grid on;
    
    subplot(3,2,4);
    plot(1:iterations, FBL_vals, 'g');
    xlabel('Iteration');
    ylabel('FBL');
    title('FBL over Iterations');
    ylim([0 1.2]);
    grid on;
    
    subplot(3,2,6);
    plot(1:iterations, FBR_vals, 'b');
    xlabel('Iteration');
    ylabel('FBR');
    title('FBR over Iterations');
    ylim([0 1.2]);
    grid on;
    
    set(gcf, "Position", [100, 100, 1200, 400]);
    
    % Add text box with parameter values
    paramText = sprintf('Parameters:\n\\alpha = %.2f\n\\beta = %.2f\n\\gamma = %.2f\n\\zeta = %.2f', ...
                        alpha, beta, gamma, zeta);
    annotation(gcf, 'textbox', [0.015, 0.85, 0.2, 0.1], ...
               'String', paramText, ...
               'FitBoxToText', 'on', ...
               'BackgroundColor', 'white', ...
               'EdgeColor', 'black', ...
               'FontSize', 10);
    
    filename = sprintf('%d_BL_BR.png', test);
    saveas(gcf, fullfile(path, filename));
end

close all;