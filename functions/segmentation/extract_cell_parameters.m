function [params] = extract_cell_parameters(I,I_seg)
%EXTRACT_CELL_PARAMETERS, Stefano Travaglino, Zhu lab, 2020
%--------------------------------------------------------------------------
%OVERVIEW: use segmentation mask to extract useful parameters from each
%cell segmented.

%INPUTS:
%I - original raw image
%I_seg = segmentation mask: 0 in background, i through n where cells are,
%where n is number of cells in the image.

%OUTPUTS:
%params: n by n_params array, where each row is [cell_id, x, y, MFI]
%--------------------------------------------------------------------------

n_cells = max(max(I_seg));

%loop through each cell (they are in order from top left bottom right)
params = zeros(n_cells,4);
for j=1:n_cells
    if j==81
        j=j;
    end
    cell_mask = I_seg==j;
    centroid = find_centroid(cell_mask);
    MFI = mean(mean(I(cell_mask)));
    %save results
    % params = [id, x, y, MFI]
    params(j,:)=[double(j), centroid, MFI];
end
end