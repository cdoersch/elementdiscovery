function pyramid = constructFeaturePyramid(data, params, levels)
% levels: What level of pyramid to compute the features for.
%
% Author: saurabh.me@gmail.com (Saurabh Singh).

if nargin < 3
  levels = [];
end
%imPaths = getImagePaths(data, imgsHome);
pyramid = constructFeaturePyramidForImg(data, params, levels);
end
