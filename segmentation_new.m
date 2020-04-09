file_name = 'C:\Users\stefano\Desktop\ca_flux_analysis\example_data\c_1.tif';
I = imread(file_name,500);
% [path,name,ext] = fileparts(app.file);

% Initialize parameters
% Set to default if the input is empty!!!!
ThreshSensitivity = 0.4;
ThreshWindowR = 32;
AreaMulti = 400;
EccentricityMulti = 0.65;
SolidityMulti = 0.5;
AreaCleanLow = 90;
AreaCleanHigh = 600;
EccentricityClean = 0.9;
SolidityClean = 0.8;            

% app.regionStats = cell(app.nFrames);
seD = cell(6,1);
for i=1:length(seD)
    seD{i} = strel('disk',i);
end

% for frameNo=1:app.nFrames
    
    %load image
%     seed = app.img(:,:,frameNo);

%generate global cell mask
T = adaptthresh(I,ThreshSensitivity,'neigh',[2*ThreshWindowR+1 2*ThreshWindowR+1]);
bw = imbinarize(I,T);
bw = imopen(bw,seD{1});    

%Local segmentation of big region based on cell edge
cc = bwconncomp(bw,8);
stats = regionprops(cc,'Area','Eccentricity','Solidity');
idx = find([stats.Area] > AreaMulti | [stats.Eccentricity] > EccentricityMulti ...
    | [stats.Solidity] < SolidityMulti); %candidates of multiple cells
if ~isempty(idx)        
    I1 = imsharpen( I,'Radius',3,'Amount',1);
    bw1 = ismember(labelmatrix(cc), idx); %subregion
    I1(~bw1) = 0;

    %detect strong local maximas
    maxima = regionmaxima(I1,seD{5});
    maxima = imopen(maxima,seD{1});
    maxima = maxima & bw1;

    %detect weak local maximas
    I2 = I1;
    for tarNo=1:5
        I2(imdilate(maxima,seD{5})) = 0;
        newmaxima = combinedmaxima(I2,bw1,maxima,seD);
        if isempty(nonzeros(newmaxima ~= maxima))        
            break
        end    
        maxima = newmaxima;    
    end

    %watershed segmentation based on local maximas
    D = -bwdist(~maxima,'chessboard');
    L = watershed(D,8);
    ridge = L == 0;
    ridge(~bw1) = 0;
    bw(ridge) = 0;
end

%clean up cell mask
cc = bwconncomp(bw,8);
stats = regionprops(cc,'Area','Eccentricity','Solidity');
idx = find([stats.Area] > AreaCleanLow & [stats.Area] < AreaCleanHigh & ...
    [stats.Eccentricity] < EccentricityClean & [stats.Solidity] > SolidityClean);
bw = ismember(labelmatrix(cc), idx);
bw = imclearborder(bw);
border = bwperim(bw);
% border = imdilate(border,seD{1});
%     app.imgSeg(:,:,frameNo,:) = imoverlay(squeeze(app.imgDisp(:,:,frameNo,1)),border,'r');
overlay = imoverlay(I,border,[.3 1 .3]);
imshow(overlay)
%Calculate cell statistics
cc = bwconncomp(bw,8);
app.regionStats{frameNo} = regionprops(cc,seed,'Area','Centroid','MeanIntensity', ...
    'Image','BoundingBox','Eccentricity','Solidity','PixelIdxList');            


close(d)

% Display segmentated images
app.SegmentationButton.Enable = 'on';
app.TrackingButton.Enable = 'off';
app.ImageSelection.SelectedObject = app.SegmentationButton;
app.FrameSlider.Value = app.nFrames;
app.Image.ImageSource = squeeze(app.imgSeg(:,:,frameNo,:));
