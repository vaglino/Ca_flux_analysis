function aligned_tracks = align_tracks(tracks,threshold)
% ALIGN_TRACKS - Stefano Travaglino, Zhu lab 2020
%--------------------------------------------------------------------------
% DESCRIPTION: align tracks by the time each cell lands on the surface

%INPUTS:
% tracks - cell array containing all tracks with parameters
% threshold - minimum time a cell has been in the field of view

%OUTPUT:
%aligned_tracks - tracks of length=threshold that are aligned and normalized
%by resting fluorescence
%--------------------------------------------------------------------------
n_tracks = numel(tracks); %total number of cells 
n_params = size(tracks{1},2); %n of columns in extracted params matrix
n_frames = size(tracks{1},1); %n of rows in params matri

plotting = true; %enable plotting in this function

if plotting
    figure
    hold on
end

aligned_tracks = cell(n_tracks,1);

for i=1:n_tracks
    
    track = tracks{i};
    start_ind = find(track(:,1) ~= 0, 1); %find first non 0 index
    %align track from first non zero index
    aligned_track = track(start_ind-1 : start_ind+threshold-1,:);  
    
    %proper normalization
    aligned_track(:,4) = normalize_fluorescence(aligned_track(:,4)); 
    
    if plotting
        plot(aligned_track(:,4),'linewidth',1)
    end
    aligned_tracks{i} = aligned_track;
    
end

tracks_matrix = cat(3, aligned_tracks{:});   %convert cell array into 3D matrix
means_matrix = mean(tracks_matrix, 3); %and get the average across 3rd dimension gg
sem_matrix = std(tracks_matrix, [], 3)./sqrt(size(tracks_matrix,3)); %get standard error

figure
plot(means_matrix(:,4))
e=errorbar(1:10:length(means_matrix(:,4)),...
            means_matrix(1:10:end,4),...
            sem_matrix(1:10:end,4),...
            'linewidth',1.5);
e.LineWidth = 1;

end

function norm_track = normalize_fluorescence(track_MFI)
%normalize track by resting fluorescence recorded just after the cell has
%landed on the surface

delay = 20; %may change this for better results

MFI_smooth = smooth(track_MFI,20);
dIdt = diff(MFI_smooth);
[~, max_i] = max(dIdt(2:50)); %skip first index, because there's a jump 
max_i = max_i+1;

I0 = MFI_smooth(max_i+delay);

norm_track = MFI_smooth./I0;

% subplot(2,2,1)
% plot(dIdt)
% subplot(2,2,2)
% plot(MFI_smooth)
% subplot(2,2,3)
% plot(norm_track)
end