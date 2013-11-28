function [pick,suppressclust] = myNms(boxes, overlap)
% top = nms_fast(boxes, overlap)
% Non-maximum suppression. (FAST VERSION)
% Greedily select high-scoring detections and skip detections
% that are significantly covered by a previously selected
% detection.
% NOTE: This is adapted from Pedro Felzenszwalb's version (nms.m),
% but an inner loop has been eliminated to significantly speed it
% up in the case of a large number of boxes
% Tomasz Maliseiwicz (tomasz@cmu.edu)

if isempty(boxes)
  pick = [];
  return;
end

x1 = boxes(:,1);
y1 = boxes(:,2);
x2 = boxes(:,3);
y2 = boxes(:,4);
s = boxes(:,end);

area = (x2-x1+1) .* (y2-y1+1);
[~, I] = sort(s);
if(nargout>=2)
  suppressclust=zeros(size(I));
end

pick = s*0;
counter = 1;
while ~isempty(I)
  
  last = length(I);
  i = I(last);  
  pick(counter) = i;
  counter = counter + 1;
  
  xx1 = max(x1(i), x1(I(1:last-1)));
  yy1 = max(y1(i), y1(I(1:last-1)));
  xx2 = min(x2(i), x2(I(1:last-1)));
  yy2 = min(y2(i), y2(I(1:last-1)));
  
  w = max(0.0, xx2-xx1+1);
  h = max(0.0, yy2-yy1+1);
  
  o = w.*h ./ area(I(1:last-1));
  f=[last; find(o>overlap)];
  if(nargout>=2)
    suppressclust(I(f))=counter-1;
  end
  
  I(f) = [];
end

pick = pick(1:(counter-1));
% top = boxes(pick,:);
