close all;

StrMap = ['#########', newline, ...
          '#+-----+#', newline, ...
          '#|#####|#', newline, ...
          '#|#+---+#', newline, ...
          '#|#|#####', newline, ...
          '#|#+---+#', newline, ...
          '#|#####|#', newline, ...
          '#+-----+#', newline, ...
          '#########'];

iterations = 20;
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

[G_loc, BV1_loc, BH1_loc, E] = parseStringMap(StrMap);
BV2_loc = BV1_loc;
BH2_loc = BH1_loc;
EV = E.*BV1_loc;
EH = E.*BH1_loc;

G_inhibition_filter = gaussianFilter2D(5, 5, 1, 1);
GV1_filter = [0.25, 0, 0; 1, 0, 0; 0.25, 0, 0];
GV2_filter = [0, 0, 0.25; 0, 0, 1; 0, 0, 0.25];
GH1_filter = [0, 0, 0; 0, 0, 0; 0.25, 1, 0.25];
GH2_filter = [0.25, 1, 0.25; 0, 0, 0; 0, 0, 0];

BV1_FB_filter = [0, 0, 0; 0, 0, 1; 0, 0, 0];
BV2_FB_filter = [0, 0, 0; 1, 0, 0; 0, 0, 0];
BH1_FB_filter = [0, 1, 0; 0, 0, 0; 0, 0, 0];
BH2_FB_filter = [0, 0, 0; 0, 0, 0; 0, 1, 0];

G = zeros(size(G_loc));
BV1_prev = zeros(size(G_loc));
BV2_prev = zeros(size(G_loc));
BH1_prev = zeros(size(G_loc));
BH2_prev = zeros(size(G_loc));
% BV1_prev = betaB*EV./(...
%     alphaB + ...
%     zetaB*EV);
% BV2_prev = betaB*EV./(...
%     alphaB + ...
%     zetaB*EV);
% 
% BH1_prev = betaB*EH./(...
%     alphaB + ...
%     zetaB*EH);
% BH2_prev = betaB*EH./(...
%     alphaB + ...
%     zetaB*EH);

BV1_history = zeros([size(BV1_prev), iterations+1]);
BV2_history = zeros([size(BV2_prev), iterations+1]);
BH1_history = zeros([size(BH1_prev), iterations+1]);
BH2_history = zeros([size(BH2_prev), iterations+1]);
G_history    = zeros([size(G), iterations+1]);

BV1_history(:,:,1) = BV1_prev;
BV2_history(:,:,1) = BV2_prev;
BH1_history(:,:,1) = BH1_prev;
BH2_history(:,:,1) = BH2_prev;
G_history(:,:,1)   = G;

for i = 1:iterations
    FBV1 = imfilter(G, BV1_FB_filter);
    FBV2 = imfilter(G, BV2_FB_filter);

    BV1 = betaB*EV.*(1 + lambdaB*FBV1)./(...
        alphaB + ...
        gammaB*BV2_prev + ...
        zetaB*EV.*(1 + lambdaB*FBV1));

    BV2 = betaB*EV.*(1 + lambdaB*FBV2)./(...
        alphaB + ...
        gammaB*BV1_prev + ...
        zetaB*EV.*(1 + lambdaB*FBV2));

    FBH1 = imfilter(G, BH1_FB_filter);
    FBH2 = imfilter(G, BH2_FB_filter);

    BH1 = betaB*EH.*(1 + lambdaB*FBH1)./(...
        alphaB + ...
        gammaB*BH2_prev + ...
        zetaB*EH.*(1 + lambdaB*FBH1));

    BH2 = betaB*EH.*(1 + lambdaB*FBH2)./(...
        alphaB + ...
        gammaB*BH1_prev + ...
        zetaB*EH.*(1 + lambdaB*FBH2));

    G_in = imfilter(BV1, GV1_filter) + ...
        imfilter(BV2, GV2_filter) + ...
        imfilter(BH1, GH1_filter) + ...
        imfilter(BH2, GH2_filter);

    G = betaG*G_in./( ...
        alphaG + ...
        gammaG*imfilter(G, G_inhibition_filter) + ...
        zetaG*G_in);

    BV1_history(:,:,i+1) = BV1;
    BV2_history(:,:,i+1) = BV2;
    BH1_history(:,:,i+1) = BH1;
    BH2_history(:,:,i+1) = BH2;
    G_history(:,:,i+1)   = G;

    BV1_prev = BV1;
    BV2_prev = BV2;
    BH1_prev = BH1;
    BH2_prev = BH2;
end

BOSV = BV1 - BV2;
BOSH = BH1 - BH2;

BOS = BOSH + BOSV;


figure; imagesc(E); axis image; colorbar;
figure; imagesc(BV1_loc); axis image;
figure; imagesc(BH1_loc); axis image;
figure; imagesc(G_loc); axis image;
figure; imagesc(BOS); axis image; colorbar;


horizontal_locations = [2,5; 8,5; 4,5; 4,6; 4,7];
vertical_locations = [5,2; 5,4];
g_locations = [3,5; 5,5; 7,5; 5,3];  % <-- Define G locations as needed

numVLocations = size(vertical_locations, 1);
numHLocations = size(horizontal_locations, 1);
numGLocations = size(g_locations, 1);
maxPerRow = 5;

numVRows = ceil(numVLocations / maxPerRow);
numHRows = ceil(numHLocations / maxPerRow);
numGRows = ceil(numGLocations / maxPerRow);
iterations_axis = 1:iterations+1;

%% --- Vertical B Cells ---
figure;
for idx = 1:numVLocations
    row = vertical_locations(idx, 1);
    col = vertical_locations(idx, 2);

    BV1_trace = squeeze(BV1_history(row, col, :));
    BV2_trace = squeeze(BV2_history(row, col, :));

    subplot(numVRows, min(maxPerRow, numVLocations), idx);
    plot(iterations_axis, BV1_trace, 'b-', 'LineWidth', 1.5); hold on;
    plot(iterations_axis, BV2_trace, 'r--', 'LineWidth', 1.5);
    title(sprintf('Vertical (%d, %d)', row, col));
    xlabel('Iter'); ylabel('Act.');
    legend('BV1','BV2','Location','southeast');
    grid on;
end
sgtitle('Vertical B Cell Activity Over Iterations');

%% --- Horizontal B Cells ---
figure;
for idx = 1:numHLocations
    row = horizontal_locations(idx, 1);
    col = horizontal_locations(idx, 2);

    BH1_trace = squeeze(BH1_history(row, col, :));
    BH2_trace = squeeze(BH2_history(row, col, :));

    subplot(numHRows, min(maxPerRow, numHLocations), idx);
    plot(iterations_axis, BH1_trace, 'b-', 'LineWidth', 1.5); hold on;
    plot(iterations_axis, BH2_trace, 'r--', 'LineWidth', 1.5);
    title(sprintf('Horizontal (%d, %d)', row, col));
    xlabel('Iter'); ylabel('Act.');
    legend('BH1','BH2','Location','southeast');
    grid on;
end
sgtitle('Horizontal B Cell Activity Over Iterations');

%% --- G Cells ---
figure;
for idx = 1:numGLocations
    row = g_locations(idx, 1);
    col = g_locations(idx, 2);

    G_trace = squeeze(G_history(row, col, :));

    subplot(numGRows, min(maxPerRow, numGLocations), idx);
    plot(iterations_axis, G_trace, 'k-', 'LineWidth', 1.5);
    title(sprintf('G Cell (%d, %d)', row, col));
    xlabel('Iter'); ylabel('Act.');
    legend('G','Location','southeast');
    grid on;
end
sgtitle('G Cell Activity Over Iterations');
