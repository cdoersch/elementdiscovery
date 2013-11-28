% greedily select elements based on coverage increase at a given level of purity.
% indata: is a set of detections in standard format, or a cell array where each
%         cell contains detections.  if a cell array, it is assumed that the
%         detection scores are not comparable, and so the purity-based selection
%         of detections is done separately for the detections in each cell.
% ispos: a logical array specifying which images are positive.  indata(:,7) should 
%        index into ispos.
% thresh: the purity threshold. For a given element, all detections are ranked by
%         score, and purity is computed for the first up through the k'th detection,
%         for all k.  We take all detections up to the last point in the ranking
%         when the purity is above thresh.
% ntosel: the number of elements to select.
% conf: a struct with the additional optional parameters:
%  - 'useoverlap': if true, compute coverage in terms of pixels
%                  (as was done for the purity-coverage plot).
%                  if false [default], compute coverage based on the number of
%                  firings on the positive set, where each firing is weighted
%                  by 1/(n), where n is the number of other detections that the
%                  patch overlaps with (jaccard overlap >.5).  Hence it's a 'softened'
%                  version of the greedy deduplication used in the SIGGRAPH 12 paper.  
%                  This is the ranking used for the indoor67 experiments, and to be 
%                  honest it doesn't work very well.  
function [res,coverageinc]=greedySelectDetrsCoverage(indata,ispos,thresh,ntosel,conf)
try
  if(~iscell(indata))
    indata={indata};
  end
  if(~exist('conf','var'))
    conf=struct();
  end
  % first select the detections above the given level of purity.
  % Do it separately for each cell of indata.
  for(m=1:numel(indata))
    [data,elid]=distributeby(indata{m},indata{m}(:,6));
    ispos=ispos(:)'>0;
    for(i=1:numel(data))
      dat=data{i};
      [~,ord]=sort(dat(:,5),'descend');
      dat=dat(ord,:);
      isposdat=ispos(dat(:,7));
      % compute purity at each point in the ranking.
      purity=cumsum(isposdat)./(1:numel(isposdat));
      tokeep=max(find(purity>=thresh));
      if(purity(end)>thresh)
        disp(['warning: min purity for ' num2str(elid(i)) ' was ' num2str(purity(end))]);
      end
      if(~dsbool(conf,'useoverlap'))
        data{i}=dat(1:tokeep,[1:4 6:7]);
        imind=6;
      else
        data{i}=dat(1:tokeep,:);
        imind=7;
      end
      data{i}=data{i}(ispos(data{i}(:,imind)),:);
      [~,ord]=sort(data{i}(:,imind));
      data{i}=data{i}(ord,:);
      if(mod(i,100)==0)
        disp(i);
      end
      if(dsbool(conf,'legaldetrs'))
        data{i}(~ismember(data{i}(:,imind-1),conf.legaldetrs),:)=[];
      end
    end
    data=cell2mat(data);
    indata{m}=data;
  end
  if(~dsbool(conf,'useoverlap'))
    % compute based on pixels
    indata=cell2mat(indata(:));
    [~,indata(:,imind)]=ismember(indata(:,imind),unique(indata(:,imind)));
    [res,coverageinc]=greedySelectDetrsCoveragemex(int64(indata),int64(ntosel));
  else
    % compute based on detection counts.
    indata=cell2mat(indata(:));
    maxdetr=max(indata(:,6));
    indata=distributeby(indata,indata(:,7));
    selected=[];
    % For each image, compute a matrix specifying which bounding boxes overlap with
    % which other bounding boxes.
    for(i=1:numel(indata))
      indata{i}=indata{i}(myNmsClass(indata{i},.5),:);
      ovlp{i}=computeOverlap(indata{i}(:,1:4),indata{i}(:,1:4),'pedro')>.5;
    end
    for(j=1:ntosel)
      % cts accumulates, for each detector, the total number of detections by that
      % detector, weighted by the number of overlaps.
      cts=zeros(maxdetr,1);
      for(i=1:numel(ovlp))
        % for each detection d in a given image, sovl(d) is the inverse of the number
        % of other detections (from other previously-selected detectors) that d overlaps with.
        sovl=1./(1+sum(ovlp{i}(ismember(indata{i}(:,6),selected),:),1));
        % Add the values in sovl to the per-detector counts in cts.  Note this
        % can't be vectorized because one image may have repeated detectors.
        for(t=1:numel(sovl))
          cts(indata{i}(t,6))=cts(indata{i}(t,6))+sovl(t);
        end
      end
      cts(selected)=-Inf;
      [coverageinc(j),selected(j)]=max(cts);
      disp(['selected ' num2str(selected(j)) ' count ' num2str(coverageinc(j))]);
    end
    res=selected(:);
  end
  %keyboard
catch ex,dsprinterr;end
end
