% Author: Carl Doersch (cdoersch at cs dot cmu dot edu)
% a relatively efficient function to convert a HOG pyramid index 
% (i.e. pyramid level and grid cell) into a position
% in pixel space.
function metadata=pyridx2pos(idx,pyrsc,patchcellsize,pyramid)%pyrcanosz,pyrsc,prSize,pcSize,sBins,imsize);
pyrcanosz=pyramid.canonicalScale;
sBins=pyramid.sbins;
imsize=[pyramid.canonicalSize.nrows pyramid.canonicalSize.ncols];
prSize=patchcellsize(1);
pcSize=patchcellsize(2);
levSc = pyramid.scales(pyrsc);
levSc=levSc(:);
canoSc = pyrcanosz;
% [x1 x2 y1 y2]
levelPatch = [ ...
  zeros(size(levSc)), ...
  round((prSize + 2) .* sBins .* levSc ./ canoSc) - 1, ...
  zeros(size(levSc)), ...
  round((pcSize + 2) .* sBins .* levSc ./ canoSc) - 1, ...
];
x1=idx(:,2);
y1=idx(:,1);
xoffset = floor((x1 - 1) .* sBins .* levSc / canoSc) + 1;
  yoffset = floor((y1 - 1) .* sBins .* levSc / canoSc) + 1;
  thisPatch = bsxfun(@plus,levelPatch, [xoffset xoffset yoffset yoffset]);

  metadata.x1 = max(1,thisPatch(:,1));
  metadata.x2 = min(thisPatch(:,2),imsize(2));
  metadata.y1 = max(1,thisPatch(:,3));
  metadata.y2 = min(thisPatch(:,4),imsize(1));
%if(metadata.x2>936)
%keyboard
%end
end

