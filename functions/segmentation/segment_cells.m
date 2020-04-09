function [I_seg] = segment_cells(I)


%adaptive contrast
I_eq = adapthisteq(I); 
% subplot(2,3,2)
% imshow(I_eq)

%median filter to get rid of the noise
I_eq = medfilt2(I_eq);
%get rid of objects touching the borders as they may be half cells or
%artifacts
I_eq = imclearborder(I_eq);
%smooth image with a kernel of 10x10 size
% windowWidth = 3;
% kernel = ones(windowWidth) / windowWidth .^ 2;
% I_smooth = imfilter(I_eq, kernel, 'replicate');

%binarize greyscale image with adaptive threshold
bw = imbinarize(I_eq, 'adaptive','ForegroundPolarity','bright');
% figure
% imshow(bw)

% binarizie after smoothing (doesn't work as well, sizes too big)
% bw1 = imbinarize(I_smooth, 'adaptive','ForegroundPolarity','bright');
% figure
% imshow(bw1)


bw2 = imfill(bw,'holes'); %fill holes
% subplot(2,3,5)
% imshow(bw2)
bw3 = imopen(bw2, ones(2,2)); %topologically open image
% subplot(2,3,6)
% imshow(bw3)
bw4 = bwareaopen(bw3, 30); %get rid of spots smaller than 30 px (empirical)
% subplot(2,3,6)
% imshow(bw4)
bw4 = imclearborder(bw4,4); 
% imshow(bw4)

%sanity check
bw4_perim = bwperim(bw4); 
overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
% figure
% imshow(overlay1)

% 
% % code for finding brightest spot inside cells (usefulness questionable)
% % figure(2)
% mask_em = imextendedmax(I_smooth, 5);
% % subplot(2,3,1)
% % imshow(mask_em)
% 
% mask_em = imclose(mask_em, ones(5,5));
% % subplot(2,3,2)
% % imshow(mask_em)
% mask_em = imfill(mask_em, 'holes');
% % subplot(2,3,3)
% % imshow(mask_em)
% mask_em = bwareaopen(mask_em, 10);
% % subplot(2,3,4)
% % imshow(mask_em)
% overlay2 = imoverlay(I_smooth, bw4_perim | mask_em, [.3 1 .3]);
% % imshow(overlay2)
% 
% I_smooth_c = imcomplement(I_smooth);
% % subplot(2,3,5)
% % imshow(I_smooth_c)

% find distance from the border of each cell
D = bwdist(~bw4);
% imshow(D,[])
D = -D;

%smooth distances (gives better result empirically)
windowWidth = 5;
kernel = ones(windowWidth) / windowWidth .^ 2;
D = imfilter(D, kernel, 'replicate');

% imshow(D,[])
% title('Complement of Distance Transform')

% watershed algorithm to segment cells
L = watershed(D);
L(~bw4) = 0;

% figure
% rgb = label2rgb(L,'jet',[.5 .5 .5]);
% imshow(rgb)
% title('Watershed Transform')

I_seg = L;
end