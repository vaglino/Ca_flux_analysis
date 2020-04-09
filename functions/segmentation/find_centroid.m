function [centroid] = find_centroid(image)

[y, x] = ndgrid(1:size(image, 1), 1:size(image, 2));
centroid = [mean(x(logical(image))), mean(y(logical(image)))];

end
