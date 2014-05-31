try
dsload('ds.batchfordetr');
dsload('ds.classperbatch');
imgs=dsload('ds.imgs{ds.conf.currimset}');
% load the detections that were used for training on the previous round.
% prevdets countains the bounding boxes; prevweights contains the alpha's
% from the paper, and prevfeats is the actual feature vectors.
prevdets=dsload(['ds.round.prevdets{' num2str(dsidx) '}'],'clear');
dsload(['ds.round.prevweights{' num2str(dsidx) '}']);
if(isfield(ds.round,'prevweights'))
  prevweights=ds.round.prevweights{dsidx};
  ds.round=rmfield(ds.round,'prevweights');
else
  prevweights=ones(size(prevdets,1),1);
end
dsload('ds.round.roundid');
% if there's a local directory, load the features from there (we're guaranteed
% to be the only machine optimizing this detector).  Otherwise they're
% stored in dswork.
if(dsfield(ds,'sys','distproc','localdir'))
  if(ds.round.roundid>5)
    % if you're running out of disk space you can delete data from old rounds...but note 
    % that this will make it impossible to back up.
    delete([ds.sys.distproc.localdir 'prevfeats' num2str(dsidx) '_' num2str(ds.round.roundid-2) '.mat']);
  end
  load([ds.sys.distproc.localdir 'prevfeats' num2str(dsidx) '_' num2str(ds.round.roundid-1) '.mat']);
else
  prevfeats=dsload(['ds.round.prevfeats{' num2str(dsidx) '}'],'clear');
end

dsload('ds.round.myiminds');
mydetrs=ds.batchfordetr(ds.batchfordetr(:,2)==dsidx,1);

% Discard any irrelevant patches from previous rounds.  That generally means
% any patch from an image that we re-ran detection on during hte current round.
% The first patch in the list for each detector, however, is generally 
% kept (control this behavior with the ds.round.exceptfirst flag).  This patch
% is the one that was randomly sampled to initialize the cluster.  This behavior
% is an artifact from versions of this codebase which treated the initial patch
% as special and broke if it disappeared.  It's probably not needed, but
% I had this flag set for the experiments in the paper.
dsload('ds.round.discardprevpatches');
dsload('ds.round.exceptfirst');
tokeep=~ismember(prevdets(:,7),ds.round.myiminds);
[~,candidatepatches]=ismember(mydetrs,prevdets(:,6),'R2012a');
if(~all(imgs.label(prevdets(candidatepatches,7))==ds.classperbatch(dsidx)))
  error('classperbatch wrong');
end
tokeep(candidatepatches)=true;
if(dsbool(ds.round,'discardprevpatches'))
  tokeep(:)=false;
  if(dsbool(ds.round,'exceptfirst'))
    tokeep(candidatepatches)=true;
  end
end
discardifnew=prevdets(candidatepatches(tokeep(candidatepatches)),[6:7]);
dets={prevdets(tokeep,:)};
feats={prevfeats(tokeep,:)};
allovlweight={prevweights(tokeep,:)};
clear prevdets;
clear prevfeats;
clear prevweights;

% Load the detections from the current round of detection.
for(i=1:numel(ds.round.myiminds))
  if(~isempty(ds.round.newfeat{dsidx,i}.assignedidx))
    tokeep=~ismember(ds.round.newfeat{dsidx,i}.assignedidx(:,6:7),discardifnew,'rows');
    dets{end+1}=ds.round.newfeat{dsidx,i}.assignedidx(tokeep,:);
    feats{end+1}=double(ds.round.newfeat{dsidx,i}.feat(tokeep,:));
    if(isfield(ds.round.newfeat{dsidx,i},'ovlweights'))
      allovlweight{end+1}=ds.round.newfeat{dsidx,i}.ovlweights(tokeep,:);
    else
      allovlweight{end+1}=ones(size(feats{end},1),1);
    end
  end
  if(mod(i,100)==0)
    disp(['img ' num2str(i) ' of ' num2str(numel(ds.round.myiminds))])
  end
end
ds.newdets{dsload('ds.round.roundid'),dsidx}=structcell2mat(dets(2:end)');

dets=structcell2mat(dets(:));
allovlweight=structcell2mat(allovlweight(:));
feats=structcell2mat(feats(:));
if(size(feats,1)>500000)
  error('featsall too big')
end

% Distribute the detections in this batch by detector id, so we can train one
% detector at a time.
[dets feats allovlweight idforcell]=distributeby(dets, feats, allovlweight, dets(:,6));
if(~all(idforcell==mydetrs(:)))
  idforcell
  mydetrs
  error('something got out of order!');
end

% load the actual detectors.
ctrs=dsload(['ds.round.detectors{' num2str(dsidx) '}'],'clear');
newctrs=zeros(size(ctrs));
resdets={};

ds.round.newfeat={};
nsv=[];
for(i=1:numel(mydetrs))
  a=tic;
  mymemory;
  weights=allovlweight{i};

  disp(['optimizing:' num2str(mydetrs(i))]); 
  disp(['total features:' num2str(size(feats{i},1))]);

  % Pull out the detector and optimize it.  See doGradientDescentproj.
  ctr=effstridx(ctrs,i);
  [newctrtmp,scores]=ds.conf.params.graddescfun(feats{i}',imgs.label(dets{i}(:,7))==ds.classperbatch(dsidx),[ctr.w ctr.b]',weights,dsload('ds.round.roundid'));
  newctr{i,1}=ctr;
  newctrtmp=newctrtmp(:)';
  newctr{i}.w=newctrtmp(1:end-1);
  newctr{i}.b=newctrtmp(end);

  dets{i}(:,5)=scores(:);

  % discard any detections that have low scores.  We keep any detections
  % with score higher than -.02/round_id, and keep at least as many detections
  % as the number of dimensions divided by 5.  Any less than this and the detections
  % on the next round could cause the element to overfit to a tiny number of patches
  % on the next round, and the round after that it will fire all over the place.
  thr=sort(scores,'descend');
  thr=min(-.02/dsload('ds.round.ndetrounds'),thr(min(ceil(size(ctr.w,2)/5),numel(thr))));
  if(scores(1)<0&&dsload('ds.round.roundid')<4),error('candidate patch died early.');end
  scores(1)=Inf;%make sure we keep the first one, since the rest of the code assumes it's there.
  feats{i}=feats{i}((scores>=thr)',:);
  dets{i}=dets{i}(scores>=thr,:);
  allovlweight{i}=allovlweight{i}(scores>=thr);
  nsv(i,1)=sum(scores>=thr);
  toc(a)
end

% convert the set of detections and features into giant matrices to save.
dets=cell2mat(dets(:));
feats=cell2mat(feats(:));
ds.nextround.prevdets{dsidx}=dets;
ds.nextround.nsv{dsidx}=nsv;
% again, if we have local storage, save it locally. Otherwise, save it to dswork.
if(dsfield(ds,'sys','distproc','localdir'))
  prevfeats=feats;
  save([ds.sys.distproc.localdir 'prevfeats' num2str(dsidx) '_' num2str(ds.round.roundid) '.mat'],'prevfeats');
else
  ds.nextround.prevfeats{dsidx}=feats;
end
ds.nextround.detectors{dsidx}=effstrcell2mat(newctr);

dssave();
ds.nextround=struct();
ds.round=struct();
ds.newdets={};
catch ex,dsprinterr;end
