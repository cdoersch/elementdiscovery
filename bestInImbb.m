% The core of the detection code.  Don't call this; call detectInIm instead.
% Note that this supports bounding boxes, but they're not used in the
% indoor67 code.
function [posall,distall,clustidall,featsall,flipall,boxidall]=bestInIm(centers,imid,conf)
    global ds;
  if(~exist('conf','var'))
    conf=struct();
  end
  conf=overrideConf(ds.conf.params,conf);
  if(~dsfield(conf,'thresh'))
    conf.thresh=-Inf;
  end
  if(numel(imid)<=2)
    imfull=im2double(getimg(imid));
  else
    imfull=im2double(imid);
  end

  noprocess=0;
  boxidall=[];
  flipall=[];
  if(dsfield(conf,'detsforclass'))% the flag indicating that we're using Pascal bounding boxes.
    annot=getannot(imid);
    bbs=[annot.x1 annot.y1 annot.x2 annot.y2];
    classes=[annot.label];
    occl=annot.occluded;
    difficult=annot.difficult;
    boxid=annot.boxid;
    flip=zeros(size(bbs,1),1);
    bbminsize=[bbs(:,4)-bbs(:,2)+1,bbs(:,3)-bbs(:,1)+1];
    if(dsbool(conf,'allowoccluded'))
      occl(:)=false;
    end
    valid=(~occl & ismember(classes,conf.detsforclass) & ~difficult & all(bsxfun(@ge,bbminsize,ds.conf.params.patchCanonicalSize),2));
    bbs(~valid,:)=[];
    boxid(~valid)=[];
    flip(~valid)=[];
  else % in this case, we simply specify one bounding box containing the full image.
    bbs=[1,1,size(imfull,2),size(imfull,1)];
    bbminsize=[bbs(:,4)-bbs(:,2)+1,bbs(:,3)-bbs(:,1)+1];
    if(~all(bbminsize>=ds.conf.params.patchCanonicalSize))
     bbs=[];
    end
    boxid=0;
    flip=0;
  end
  if(dsbool(conf,'flipall')) %add a flipped bounding box.
    bbs=[bbs;bbs];
    boxid=[boxid;boxid];
    flip=[flip;ones(size(flip))];
  end
  if(isempty(bbs))
    posall=[];
    distall=[];
    clustidall=[];
    featsall=[];
    boxid=[];
    return
  end
  for(bbidx=1:size(bbs,1))
    im=imfull(bbs(bbidx,2):bbs(bbidx,4),bbs(bbidx,1):bbs(bbidx,3),:);
    if(flip(bbidx)),im=im(:,end:-1:1,:);end
    pyramid = constructFeaturePyramid(im, ds.conf.params); % HOG feature pyramid
    pcs=round(ds.conf.params.patchCanonicalSize/ds.conf.params.sBins)-2;
    pcs(3)=size(pyramid.features{1},3);
    pcs(4)=0;
    if(numel(imid)>2)
      if(~isfield(conf,'imid'))      
        conf.imid=0
      end
    else
      conf.imid=imid;
    end
    % unentangling the fieature pyramid gives us 'features', where each row
    % is the feature vector for a single patch.  levels and indexes specify
    % where those patches were in the pyramid, and gradsums tells us the strength
    % of the gradient in each patch so we can get rid of empty ones that are likely
    % to make our detectors misfire thanks to HOG's normalization.
    [features, levels, indexes,gradsums] = unentanglePyramid(pyramid, ...
      pcs,conf);


    invalid=(gradsums<9);
    features(invalid,:)=[];
    levels(invalid)=[];
    indexes(invalid,:)=[];
    gradsums(invalid)=[];
    
    if(dsbool(conf,'multperim'))
      % findmatches finds all detections above a certain threshold
      [assignedidx, dist, clustid]=findmatches(centers,features,conf.thresh,conf);
    else
      % assigntoclosest gives us the top detection for each detector.
      [assignedidx, dist]=assigntoclosest(centers,features,1);
      if(isempty(dist))
        clustid=[];
      else
        clustid=(1:size(centers,1))';
        valid=dist>conf.thresh;
        assignedidx=assignedidx(valid);
        dist=dist(valid);
        clustid=clustid(valid);
      end
    end
    patsz=ds.conf.params.patchCanonicalSize;%allsz(resinds(k),:);
    fsz=(patsz-2*ds.conf.params.sBins)/ds.conf.params.sBins;
    imgs=getimgs();
    % convert the pyramid indexes into bounding boxes with pixel coordinates.
    pos=pyridx2pos(indexes(assignedidx,:),reshape(levels(assignedidx),[],1),fsz,pyramid);
    % if we allow more than one per image, we need to do non-maximum suppresion.
    if(dsbool(conf,'multperim'))
      pos=[pos.x1 pos.y1 pos.x2 pos.y2];
      [pos,assignedidx,dist,clustidl,clustid]=distributeby(pos,assignedidx,dist,clustid,clustid);
      for(i=1:numel(pos))
        [posinds]=myNms([pos{i} dist{i}],ds.conf.params.nmsOverlapThreshold);
        assignedidx{i}=assignedidx{i}(posinds);
        dist{i}=dist{i}(posinds);
        pos{i}=pos{i}(posinds,:);
        clustidl{i}=clustidl{i}(posinds);
      end
      assignedidx=cell2mat(assignedidx);
      dist=cell2mat(dist);
      p=cell2mat(pos);
      if(isempty(p))
        p=zeros(0,4);
      end
      clear pos;
      pos.x1=p(:,1);pos.x2=p(:,3);pos.y1=p(:,2);pos.y2=p(:,4);
      clustid=cell2mat(clustidl);
    end
    feats=features(assignedidx,:);
    pos.x1=pos.x1+bbs(bbidx,1)-1;
    pos.x2=pos.x2+bbs(bbidx,1)-1;
    pos.y1=pos.y1+bbs(bbidx,2)-1;
    pos.y2=pos.y2+bbs(bbidx,2)-1;
    if(flip(bbidx))
      medval=(bbs(bbidx,3)+bbs(bbidx,1))/2;
      tmp=medval+(medval-pos.x2);
      pos.x2=medval+(medval-pos.x1);
      pos.x1=tmp;
    end

    posall{bbidx,1}=effstr2str(pos);
    featsall{bbidx,1}=feats;
    distall{bbidx,1}=dist;
    clustidall{bbidx,1}=clustid;
    boxidall{bbidx,1}=repmat(boxid(bbidx),size(dist,1),1);
    flipall{bbidx,1}=repmat(flip(bbidx),size(dist,1),1);
  end
  posall=str2effstr(cell2mat(posall));
  featsall=cell2mat(featsall);
  distall=cell2mat(distall);
  clustidall=cell2mat(clustidall);
  boxidall=cell2mat(boxidall);
  flipall=cell2mat(flipall);
end
