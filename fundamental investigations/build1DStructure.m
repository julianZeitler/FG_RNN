function [B_in, B_out, G, W_in, W_out] = build1DStructure(StrMap, G_radius, showPlot, savePath)
    % build1DStructure - Builds 1D structure from string map
    % 
    % Inputs:
    %   StrMap   - String representation of the structure
    %   showPlot - Boolean flag to show plot (default: false)
    %   savePath - Optional path to save the figure (default: empty, no save)
    %
    % Outputs:
    %   B_in, B_out, G, W_in, W_out - Structure components
    
    if nargin < 3
        showPlot = false;
    end
    if nargin < 4
        savePath = '';
    end

    [G_loc, Bv, Bh, E] = parseStringMap(StrMap);
    B = Bv | Bh;

    % Find boundary B loop
    B_bin = B > 0;
    [B_boundaries, Map] = bwboundaries(B_bin, 8);
    B_idcs = [];
    for k = 1:length(B_boundaries)
        boundary = B_boundaries{k};
        boundary = unique(boundary, "stable", "rows");
        B_boundaries{k} = boundary;
        binary_map = zeros(size(B_bin, 1), size(B_bin, 2));
        idx = sub2ind(size(binary_map), boundary(:,1), boundary(:,2));
        binary_map(idx) = 1;
        if any(binary_map & B_bin, "all")
            B_idcs = [B_idcs, k];
        end
    end
    
    % Store B loops and create index maps
    B_loops = cell(length(B_idcs), 1);
    B_index_maps = cell(length(B_idcs), 1);
    
    % Initialize inward/outward orientation matrices
    Bv_in = zeros(size(Bv));
    Bh_in = zeros(size(Bh));
    
    % Decide for every B cell whether it's pointing inwards or outwards
    for b_idx = 1:length(B_idcs)
        loop = B_boundaries{B_idcs(b_idx)};
        n_B = size(loop, 1);
        loop_lin = sub2ind(size(B), loop(:,1), loop(:,2));
        B_loops{b_idx} = loop_lin;
        B_index_maps{b_idx} = containers.Map(loop_lin, 1:n_B);
        
        % Determine if boundary is closed (first and last points are adjacent)
        is_closed = (norm(loop(1,:) - loop(end,:)) <= sqrt(2)+0.0001);
        
        % Determine boundary orientation and inward direction
        if is_closed
            for i = 1:n_B
                current_pos = loop(i, :);
                y = current_pos(1);
                x = current_pos(2);
                
                if Bv(y, x) == 1
                    if Map(y, x+1) > Map(y,x)
                        Bv_in(y,x) = 1;
                    else
                        Bv_in(y,x) = -1;
                    end
                end
                if Bh(y, x) == 1
                    if Map(y-1, x) > Map(y,x)
                        Bh_in(y,x) = 1;
                    else
                        Bh_in(y,x) = -1;
                    end
                end
            end
        else
            Bh_in = Bh;
            Bv_in = Bv;
            % y_1 = loop(1,1);
            % x_1 = loop(1,2);
            % if Bv(y_1, x_1) == 1
            %     Bv_in(y_1, x_1) = 1;
            % end
            % if Bh(y_1, x_1) == 1
            %     Bh_in(y_1, x_1) = 1;
            % end
            % 
            % for i = 2:n_B
            %     prev_pos = loop(i-1, :);
            %     current_pos = loop(i, :);
            %     y = current_pos(1);
            %     x = current_pos(2);
            % 
            %     direction = current_pos - prev_pos;
            % 
            %     if direction(1) < 0 && direction(2) < 0
            % 
            %     elseif direction(1) < 0 && direction(2) == 0
            %         if Bv(y, x) == 1 && Bh(y, x) == 1
            %         elseif Bv(y, x) == 1
            %         elseif Bh(y, x) == 1
            %         end
            %     elseif direction(1) < 0 && direction(2) > 0
            %     elseif direction(1) == 0 && direction(2) < 0
            %     elseif direction(1) == 0 && direction(2) > 0
            %     elseif direction(1) > 0 && direction(2) < 0
            %     elseif direction(1) > 0 && direction(2) == 0
            %     elseif direction(1) > 0 && direction(2) > 0
            %     end
            % end
        end
    end

    [G_y, G_x] = find(G_loc);
    G_y = G_y(:);
    G_x = G_x(:);
    G_yx = [G_y, G_x];
    n_total_G = length(G_y);

    G_indices = zeros(size(B));
    idx = 1;
    for y_idx = 1:size(B, 1)
        for x_idx = 1:size(B, 2)
            if ismember([y_idx, x_idx], G_yx, 'rows')
                G_indices(y_idx, x_idx) = idx;
                idx = idx + 1;
            end
        end
    end

    % Initialize outputs
    B_in = cell(length(B_loops), 1);
    B_out = cell(length(B_loops), 1);
    W_in = cell(length(B_loops), 1);
    W_out = cell(length(B_loops), 1);
    
    for shape_idx = 1:length(B_loops)
        n_B = length(B_loops{shape_idx});
        B_in{shape_idx} = zeros(n_B, 1);
        B_out{shape_idx} = zeros(n_B, 1);
        W_in{shape_idx}  = sparse(n_total_G, n_B);
        W_out{shape_idx} = sparse(n_total_G, n_B);
    end
    G = zeros(n_total_G, 1);

    % Define B-cell filters (from 2D case)
    GRFs = makeGRF(G_radius, [0, pi/2]);
    filters = {GRFs{1, 1}, GRFs{1, 2}, GRFs{1, 3}, GRFs{1, 4}};
    
    for shape_idx = 1:length(B_loops)
        for g_idx = 1:n_total_G
            gy = G_yx(g_idx, 1);
            gx = G_yx(g_idx, 2);
            for f = 1:length(filters)
                F = filters{f};
                ori = (f-1)*pi/2;
                for fy = 1:size(F, 2)
                    for fx = 1:size(F, 1)
                        w = F(fy, fx);
                        if w == 0, continue; end
                        ny = gy + (fy - ceil(size(F, 1)/2));
                        nx = gx + (fx - ceil(size(F, 2)/2));
                        if ny < 1 || ny > size(B,1) || nx < 1 || nx > size(B,2)
                            continue;
                        end
                        lin_idx = sub2ind(size(B), ny, nx);
                        if ismember(lin_idx, B_loops{shape_idx})
                            g_lin_idx = G_indices(gy, gx);
                            B_index_map = B_index_maps{shape_idx};
                            b_idx = B_index_map(lin_idx);
                            if Bv(ny, nx) == 1 && Bh(ny, nx) == 1
                                if ori == 0  % Vertical orientation
                                    if Bv_in(ny, nx) == -1
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    end
                                elseif ori == pi/2  % Horizontal orientation
                                    if Bh_in(ny, nx) == -1
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    end
                                elseif ori == pi  % Vertical orientation
                                    if Bv_in(ny, nx) == 1  % Opposite direction for pi orientation
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    end
                                elseif ori == 3*pi/2  % Horizontal orientation
                                    if Bh_in(ny, nx) == 1  % Opposite direction for 3*pi/2 orientation
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w/2;
                                    end
                                end
                            elseif Bh(ny, nx) == 1  % Pure horizontal boundary
                                if ori == pi/2
                                    if Bh_in(ny, nx) == -1
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w;
                                    end
                                elseif ori == 3*pi/2
                                    if Bh_in(ny, nx) == 1  % Opposite direction for 3*pi/2 orientation
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w;
                                    end
                                end
                            elseif Bv(ny, nx) == 1  % Pure vertical boundary
                                if ori == 0
                                    if Bv_in(ny, nx) == -1
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w;
                                    end
                                elseif ori == pi
                                    if Bv_in(ny, nx) == 1  % Opposite direction for pi orientation
                                        W_in{shape_idx}(g_lin_idx, b_idx) = W_in{shape_idx}(g_lin_idx, b_idx) + w;
                                    else
                                        W_out{shape_idx}(g_lin_idx, b_idx) = W_out{shape_idx}(g_lin_idx, b_idx) + w;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    if showPlot
        % Visualization of G and B cell indices on original 2D map
        figure('Name', 'Structure Map');
        imagesc(2*E + G_loc); 
        colormap(gray); 
        axis equal tight;
        hold on;
    
        % Plot G cell indices
        for y = 1:size(G_indices, 1)
            for x = 1:size(G_indices, 2)
                g_idx = G_indices(y, x);
                if g_idx > 0
                    text(x, y, num2str(g_idx), 'Color', 'g', 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end
            end
        end
    
        % Plot B cell indices for all loops
        colors = {'r', 'b', 'm', 'c', 'y'};  % Different colors for different loops
        for loop_idx = 1:length(B_loops)
            current_loop = B_loops{loop_idx};
            color = colors{mod(loop_idx-1, length(colors)) + 1};
            
            % Convert linear indices back to subscripts
            [loop_y, loop_x] = ind2sub(size(B), current_loop);
            
            for i = 1:length(current_loop)
                y = loop_y(i);
                x = loop_x(i);
                text(x, y, sprintf('%d', i), 'Color', color, 'FontSize', 5, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        end
    
        title('G (green) and B (colored by loop) Cell Indices');
        
        % Save figure if path is provided
        if ~isempty(savePath)
            % Ensure directory exists
            [pathDir, ~, ~] = fileparts(savePath);
            if ~isempty(pathDir) && ~exist(pathDir, 'dir')
                mkdir(pathDir);
            end
                        
            % Also save as .png if the path doesn't have .fig extension
            [~, name, ext] = fileparts(savePath);
            if ~strcmp(ext, '.fig')
                % User provided a different extension or no extension
                if isempty(ext)
                    % No extension provided, save both .fig and .png
                    savefig([savePath, '.fig']);
                    print([savePath, '.png'], '-dpng', '-r300');
                else
                    % Save with the provided extension
                    print(savePath, '-dpng', '-r300');
                end
            else
                % Also create a .png version
                pngPath = fullfile(pathDir, [name, '.png']);
                print(pngPath, '-dpng', '-r300');
            end
        end
    end
end