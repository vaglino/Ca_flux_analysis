%FLUX_ANALYSIS - Stefano Travaglino, Zhu Lab, 2020
%--------------------------------------------------------------------------
%OVERVIEW: master script to load tiff stack, segment each image, track each
%cell and perform single cell flux analysis

%INPUTS:
%threshold - minimum time a cell needs to be image for, such that data is
%meaningful and useful
%--------------------------------------------------------------------------

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% load image stack
close all
%put your image pathway here
file_name = 'C:\Users\stefano\Desktop\ca_flux_analysis\example_data\c_1.tif';
tiff_info = imfinfo(file_name); % return tiff structure, one element per image

%% hyperparameters
n_frames = size(tiff_info, 1); % total number of recorded frames
n_frames = 300; 
threshold = 200; % minimum time for cells to be recorded

%% segmenting
tic
save_params = cell(n_frames,1);
points = cell(n_frames,1);
% progress_bar = waitbar(0,'Segmentation start...');
for i=1:n_frames
% for i=900:n_frames
%     if i==900
%         i=i
%     end
    if mod(i,10)==0
%         message = sprintf('Segmenting frame %d/%d',[i,n_frames]);
%         waitbar(i/n_frames,progress_bar,message);
        fprintf('Segmenting frame %d/%d\n',[i,n_frames])
    end

    %loop through all images
%     I = imread(files1(i).name);
    I = imread(file_name,i);

    %for each image segment cells
    [I_seg] = segment_cells(I);
    
    %calculate centroid and MFI for each cell
    [params] = extract_cell_parameters(I,I_seg);
    centroids = params(:,2:3);

    save_params{i} = params;
    points{i} = centroids;
end
toc
% close(progress_bar)

%% tracking

max_linking_distance = 5; %max velocity cells can have, may change this, try out what works best
max_gap_closing = 10; %number of frames where a cell is not segmented

[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing',      max_gap_closing, ...
    'Debug',              false);

parsed_tracks = parse_tracks(save_params,tracks,threshold);

%% post processing

[save, n_tot_cells] = pick_tracks(parsed_tracks,threshold);

aligned_tracks = align_tracks(save,threshold);

