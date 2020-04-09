function parsed_tracks = parse_tracks(raw_params,tracks,threshold)
%PARSE_TRACKS - Stefano Travaglino, Zhu Lab, 2020
%--------------------------------------------------------------------------
%OVERVIEW: cell id numbers after segmentation are random, therefore we need
%to use the information provided by tracking to reorganize the parameters
%obtained from segmentation. Each track provided in the input is a list of
%cell ids of a single T-cell across frames. This function uses the ids to
%extract the segmented parameters for each T-cells across frames

%INPUTS:
%raw_params - n by 1 cell array, where n is number of frames, with each array cell
%containing extracted parmeters for each cell present in that frame.
%i.e. [cell_id, x_centroid, y_centroid, MFI]

%tracks - n_tracks by 1 cell array, where n_tracks is the total number of
%T cells imaged, and each cell in the array contains a n_frame by 1 list of
%cell ids, which indicate the segmentation id of a single cell across
%frames.

%OUTPUT:
%parsed_tracks

%--------------------------------------------------------------------------
n_tracks = size(tracks,1); %total number of cells 
n_params = size(raw_params{1},2); %n of columns in extracted params matrix
n_frames = numel(raw_params);

plotting = true; %enable plotting in this function

%parse tracks, following a each single cell through its path and extracting
%its parameters from input mat

if plotting
    figure
    hold on
end

parsed_tracks = cell(n_tracks,1);

counter = 0;
for i=1:n_tracks %loop through each track
    
    track = tracks{i};
    
    %skip track if they are just NaN
    if sum(~isnan(track(:,1)))<2
        continue
    end
    
    track_params = zeros(n_frames,n_params);
    for j=1:n_frames %loop through each frame
        cell_id = track(j); %cell id in this frame
        
        %extract parameters for cell in this frame, or assign NaN if there
        %is no cell
        if isnan(cell_id)
            params = NaN(1,n_params);
        else
            params = raw_params{j}(cell_id,:); %[cell_id, x_centroid, y_centroid, MFI]
        end
        
        track_params(j,:) = params;
    end
        
    %clean track (treat NaNs at start and end, interpolate other NaNs)
    cleaned_track = clean_track(track_params,n_frames,n_params);
    
    %if plotting is enabled, plot track if it reaches minimum length
    if plotting
        if sum(cleaned_track(:,4)~=0)>=threshold && cleaned_track(1,4)==0
            plot(smooth(cleaned_track(:,4),20))
        end
    end
    
    parsed_tracks{i} = cleaned_track;
end

index_empty = cellfun(@isempty, parsed_tracks) == 0;
parsed_tracks = parsed_tracks(index_empty);

end


function [cleaned_track] = clean_track(track_params,n_frames,n_params)
% CLEAN_TRACKS, Stefano Travaglino, Zhu Lab, 2020
%--------------------------------------------------------------------------
%OVERVIEW: this function gets rid of all NaNs obtained from tracking. 
%namely it replaces NaNs at the start and at the end by 0s, and
%interpolates NaNs in the middle of the track, such that for tracks that
%were linked and had no parameter values obtained from segmentation, NaNs
%are replaced by linearly interpolated values.

%INPUTS: 
%track_params - n_frames by n_params array containing parameters for
%a single T-cell. 

%n_frames and n_params - number of frames and recorded parameters
%
%OUTPUTS:
%cleaned track - n_frames by n_params array contatining parameters for a
%single T-cell, with all NaNs replaced by appropriate values.
%--------------------------------------------------------------------------
cell_ids = track_params(:,1);
cleaned_track = track_params;

%zero pad NaNs at the beginning of each track
idx_min = find(~isnan(cell_ids), 1); %finds first non-NaN id index 
if idx_min ~= 1
    cleaned_track(1:idx_min-1,:) = 0; %sets all params to zero until cell touches glass
end

%zero pad NaNs at the end of each track
idx_max = find(~isnan(cell_ids), 1, 'last');
if idx_max ~= n_frames
    cleaned_track(idx_max+1:end,:) = 0;
end

%interpolate NaNs inside track
nan_idx = isnan(cleaned_track(:,1)); %detremine indexes of NaN frames
frames= (1:n_frames)';
nnan_frames = frames(~nan_idx); %indexes of non NaN frames
nan_frames = frames(nan_idx); %indexes of NaN frames

for i=2:n_params %interpolate values for all parameters other than cell id
    parameter_vals = track_params(~nan_idx,i);
    cleaned_track(nan_idx,i) = interp1(nnan_frames, parameter_vals, nan_frames);
end


end