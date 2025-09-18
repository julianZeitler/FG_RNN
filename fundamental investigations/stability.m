beta_B = 1;
zeta_B = 0.9;
gamma_SE = 1;
B_max = beta_B/zeta_B;

g_relu      = @(B) max(0, B);
g_relu_d    = @(B) double(B > 0);

g_shiftrelu   = @(B) max(0, B - 0.5);
g_shiftrelu_d = @(B) double(B > 0.5);

g_sigmoid   = @(B) 1 ./ (1 + exp(-(B - 2)));
g_sigmoid_d = @(B) g_sigmoid(B) .* (1 - g_sigmoid(B));

E = 1;

B = linspace(0, B_max*1.2, 1000);

G = {g_relu, g_shiftrelu, g_sigmoid};
Gd = {g_relu_d, g_shiftrelu_d, g_sigmoid_d};
titles = {'ReLU', 'Shifted ReLU', 'Sigmoid'};

figure;
for i = 1:3
    g = G{i};
    g_prime = Gd{i};

    f1 = zeta_B * (E + gamma_SE * g(B));
    f2 = (beta_B + zeta_B * B) .* gamma_SE .* g_prime(B);

    subplot(1, 3, i);
    plot(B, f1, 'b-', 'LineWidth', 2); hold on;
    plot(B, f2, 'r--', 'LineWidth', 2);
    xlabel('B'); ylabel('Value');
    title(['g(B): ', titles{i}]);
    legend('f_1(B)', 'f_2(B)', 'Location', 'Best');
    grid on;
end

set(gcf, "Position", [100, 100, 1800, 400]);
    
% Add text box with parameter values
paramText = sprintf("Parameters:\n\\beta_B = %.2f\n\\gamma_{SE} = %.2f\n\\zeta_B = %.2f\n\nFunctions:\nf_1(B) = \\zeta_B(E + \\gamma_{SE}g(B))\nf_2(B) = (\\beta_B + \\zeta_BB)\\gamma_{SE}g'(B)", ...
                    beta_B, gamma_SE, zeta_B);
annotation(gcf, 'textbox', [0.01, 0.87, 0.2, 0.1], ...
           'String', paramText, ...
           'FitBoxToText', 'on', ...
           'BackgroundColor', 'white', ...
           'EdgeColor', 'black', ...
           'FontSize', 10);

sgtitle('Comparison of f_1(B) and f_2(B) for Different g(B). f_1(B)>f_2(B) for stability.');