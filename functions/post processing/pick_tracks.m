function [save, n_tot_cells] = pick_tracks(raw_tracks,threshold)

% PICK_TRACKS - Stefano Travaglino, Zhu lab 3/2/2020
%--------------------------------------------------------------------------
% DESCRIPTION: gets rid of cells that were already present in the first 
% frame of recording (we don't know exactly when they landed on the 
% coverslip because they were already in the field of view at t<=0)
% then only selects tracks for cells that have been imaged for longer than
% a time threshold specified in the input

%INPUTS:
% raw_tracks - cell array containing all cleaned up tracks
% threshold - minimum time a cell has been in the field of view

%OUTPUT:
% save - cell array conatining all valid tracks of length > threshold 
%--------------------------------------------------------------------------
n_tracks = numel(raw_tracks); %total number of cells 
n_params = size(raw_tracks{1},2); %n of columns in extracted params matrix
n_frames = size(raw_tracks{1},1); %n of rows in params matri

plotting = true; %enable plotting in this function

if plotting
    figure
    hold on
end

pick_flags = ones(n_tracks,1); % 1 for acceptable, 0 not acceptable track
counter = 0;
for i=1:n_tracks
    % check if track is longer than thrshold length, as well as if cell
    % landed before t=0
    track = raw_tracks{i};
    if sum(track(:,1)~=0) <= threshold || track(1,1) ~= 0
        
        pick_flags(i) = 0; % bad track
        
    else      
        plot(smooth(track(:,4)),'linewidth',1) 
    end
end

save = raw_tracks(logical(pick_flags));
n_tot_cells = sum(pick_flags);

fprintf('%d cells were recorded in this experiment',n_tot_cells)

end
