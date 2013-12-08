% Generate a heatmap display of the classifier, which shows the pixels in
% the image that were most influential.  Basically finds which dimensions
% of the image descriptor were most influential and back-projects them to
% the max firing within each region of the spatial pyramid.  Note that this
% tends to generate a sparse representation of textured regions.  This is
% different from the visualization on the website, (that one does more to smooth
% the detections and counts all detections, rather than just the max ones).
%
% Note that you may need to adjust the scaling of the contributions on
% lines 108-109 before the heatmaps look good.  This function does not do
% that adjustment automatically because we want the heatmap values to be
% comparable across images.
%
% detr: element detectors
% im: image
% svm: final SVM
% feattransf: a function which converts max-pooled responses into a feature
% vector suitable for svm classification.
% dets: topn detections for each detector in detr
% dispimset: index of the image set where dets came from (note: in general,
%            im was loaded from the test images, dets are from training images.
% displayname: prefix for saving the displays
function [posheatmap, negheatmap] = dispClassifier(detr, im, svm, feattransf,dets,displayname,dispimset)
try

global ds;
  % first, re-generate the image descriptor.  As we find detections, add them to
  % the heat map.
    pyramid = constructFeaturePyramid(im, ds.conf.params);
    [features, levels, indexes,gradsums] = unentanglePyramid(pyramid, ...
      ds.conf.params.patchCanonicalSize/ds.conf.params.sBins-2);
  invalid=(gradsums<9);
  size(features)
  features(invalid,:)=[];
  levels(invalid)=[];
  indexes(invalid,:)=[];
  gradsums(invalid)=[];
  disp(['threw out ' num2str(sum(invalid)) ' patches']);
   patsz=ds.conf.params.patchCanonicalSize;
   fsz=(patsz-2*ds.conf.params.sBins)/ds.conf.params.sBins;
   pos=pyridx2pos(indexes,reshape(levels,[],1),fsz,pyramid);
    posy=(pos.y1 + pos.y2)/2+.000001;
    posx=(pos.x1 + pos.x2)/2+.000001;
    posheatmap=zeros(size(im(:,:,1)));
    negheatmap=zeros(size(posheatmap));
    idx=1;
    feature=zeros(5,size(detr.w,1));
    wmat=reshape(svm.w,5,[]);
    for(i=[-1 1])
      for(j=[-1 1]);
        posidx=find(i*(posy-size(im,1)/2) > 0 & j*(posx-size(im,2)/2) > 0);
        [assignedidx,dist]=assigntoclosest(detr.w,features(i*(posy-size(im,1)/2) > 0 & j*(posx-size(im,2)/2) > 0,:),1);
        dist=dist(:)-detr.b;
        posidx2=posidx(assignedidx);
        posmat=[pos.x1 pos.y1 pos.x2 pos.y2];
        posmat2=posmat(i*(posy-size(im,1)/2) > 0 & j*(posx-size(im,2)/2) > 0,:);
        if(exist('feattransf','var'))
          dist=feattransf(dist')';
        end
        feature(idx,:)=dist(:)';
        assignedidxall{idx}=assignedidx;
        posall{idx}=posmat2(assignedidx,:);
        wt=dist(:)'.*wmat(idx,:);
        posheatmap=posheatmap+genheatmap(wt(wt>0)',posmat2(assignedidx(wt>0),:),size(posheatmap));
        negheatmap=negheatmap+genheatmap(-wt(wt<0)',posmat2(assignedidx(wt<0),:),size(posheatmap));
        idx=idx+1;
      end
    end
    [feature(end,:),maxpos]=max(feature(1:end-1,:),[],1);
    wt=feature(end,:).*wmat(end,:);
    for(i=1:4)
       posheatmap=posheatmap+genheatmap(wt(wt>0&maxpos==i)',posall{i}(wt>0&maxpos==i,:),size(posheatmap));
       negheatmap=negheatmap+genheatmap(-wt(wt<0&maxpos==i)',posall{i}(wt<0&maxpos==i,:),size(posheatmap));
    end
    hm=posheatmap-negheatmap-svm.rho/numel(posheatmap);
    posheatmap=hm.*(hm>0);
    negheatmap=-hm.*(hm<0);
    % if we got detections as an argument, generate two displays to show the ones
    % that contributed most to this detection: one for negative contributions,
    % one for positive.
    if(nargin>4)
      currimset=ds.conf.currimset;
      if(exist('dispimset','var'))
        ds.conf.currimset=dispimset;%assume imgs is loaded; don't want to run dsup.
      end
      contrib=sum(feature.*wmat,1);
      [wt,todisp]=maxk(contrib,25);
      todisp(wt<0)=[];
      wt(wt<0)=[];
      todisp=detr.id(todisp);
      dsup([displayname '_pos.patchimg'],extractpatches(dets(ismember(dets(:,6),todisp),:)));
      conf=struct('dets',dets(ismember(dets(:,6),todisp),:),'detrord',todisp,...
                  'message',{cellfun(@(x) ['contribution:' num2str(x)],num2cell(wt),'UniformOutput',false)});
      mhprender('patchdisplay.mhp',[displayname '_pos.displayhtml'],conf);
      [wt,todisp]=mink(contrib,50);
      todisp(wt>0)=[];
      wt(wt>0)=[];
      todisp=detr.id(todisp);
      dsup([displayname '_neg.patchimg'],extractpatches(dets(ismember(dets(:,6),todisp),:)));
      conf=struct('dets',dets(ismember(dets(:,6),todisp),:),'detrord',todisp,...
                  'message',{cellfun(@(x) ['contribution:' num2str(x)],num2cell(wt),'UniformOutput',false)});
      mhprender('patchdisplay.mhp',[displayname '_neg.displayhtml'],conf);
      if(exist('dispimset','var'))
        ds.conf.currimset=currimset;
      end
    end
    
    feature=feature(:);
    im=repmat(rgb2gray(im),[1,1,3]);
    if(dsbool(ds.conf.params,'ovlweight'))
      weight=1000;
    else
      weight=2000;
    end
    posheatmap=heatmap2jet(posheatmap*weight)*.5+im*.5;
    negheatmap=heatmap2jet(negheatmap*weight)*.5+im*.5;
catch ex,dsprinterr;end
end

function res=heatmap2jet(heatmap)
  cmp=colormap('jet');
  res=zeros([size(heatmap) 3]);
  heatmap=round(heatmap*size(cmp,1));
  heatmap(heatmap<1)=1;
  heatmap(heatmap>size(cmp,1))=size(cmp,1);
  for(chan=1:3)
    res(:,:,chan)=reshape(cmp(heatmap,chan),size(heatmap));
  end
end
