function [G, B_vertical, B_horizontal, E] = parseStringMap(StrMap)
% parseStringMap Parses a grid string and returns binary maps for grid cells, vertical B cells and horizontal B cells.
%
%   [G, B_vertical, B_horizontal, E] = PARSEGRIDSTRING(gridString)
%
%   This function takes a multi-line string representing a grid and generates
%   three binary matrices (maps) based on specific characters:
%   1. 'G': Contains 1s where '#' characters are found.
%   2. 'B_vertical': Contains 1s where '|' or '+' characters are found.
%   3. 'B_horizontal': Contains 1s where '-' or '+' characters are found.
%
%   Input:
%     gridString - A string where each line represents a row of the grid.
%                  Lines are separated by newline characters (e.g., char(10) or '\n').
%                  Example: '#########\n#+-----+#\n...'
%
%   Outputs:
%     G     - A 2D logical (or double) array where elements are 1
%                    if the corresponding character in the gridString was '#',
%                    and 0 otherwise.
%     B_vertical - A 2D logical (or double) array where elements are 1
%                    if the corresponding character was '|' or '+', and 0 otherwise.
%     B_horizontal - A 2D logical (or double) array where elements are 1
%                    if the corresponding character was '-' or '+', and 0 otherwise.

    % Split the input string into individual lines based on newline characters.
    % strsplit returns a cell array of strings.
    lines = strsplit(StrMap, '\n');

    % Remove any empty lines that might occur, for example, if the input
    % string ends with a newline.
    lines = lines(~cellfun('isempty', lines));

    % Determine the number of rows in the grid.
    numRows = length(lines);

    % Handle empty input string case.
    if numRows == 0
        G = [];
        B_vertical = [];
        B_horizontal = [];
        return;
    end

    % Determine the number of columns. We assume all lines have the same length
    % for a rectangular grid, based on the example provided.
    numCols = length(lines{1});

    % Convert the cell array of strings into a 2D character matrix.
    % This is efficient for character-by-character operations.
    gridMatrix = char(lines);

    % Initialize the three output maps with zeros.
    % These will be populated with 1s where the conditions are met.
    G = zeros(numRows, numCols);
    B_vertical = zeros(numRows, numCols);
    B_horizontal = zeros(numRows, numCols);

    % Populate the 'G': set 1s where the character is '#'.
    G(gridMatrix == '#') = 1;

    % Populate the 'B_vertical': set 1s where the character is '|' or '+'.
    B_vertical(gridMatrix == '|' | gridMatrix == '+') = 1;

    % Populate the 'B_horizontal': set 1s where the character is '-' or '+'.
    B_horizontal(gridMatrix == '-' | gridMatrix == '+') = 1;

    E = B_vertical | B_horizontal;

end