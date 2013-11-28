% Carl Doersch (cdoersch at cs dot cmu dot edu)
% Initialize each detector in one batch.  Initial weight vectors are in
% ds.detectors{dsidx}, and initial patch features are in ds.initFeats.
% For the most part, it's just copying the data in ds.initFeats into the
% detector structure, but it also needs to set the bias.  It's important
% to give a reasonable initial estimate of the bias because if affects the
% number of patches mined in the early rounds, which can have a big impact
% on speed.

posPats=dsload('ds.initPatches');
dsload('ds.initFeats');

% Find the label for this batch and select a set of negative patches from 
% ds.initFeats. These negatives are needed to calibrate each detector in the batch.
dsload('ds.imgs');
dsload('ds.classperbatch');
dsload('ds.batchfordetr');
myclusts=ds.batchfordetr(ds.batchfordetr(:,2)==dsidx,1);
neginds=find(idxwithdefault(ds.imgs{ds.conf.currimset}.label,ds.initPatches(:,7),0)~=ds.classperbatch(dsidx));
rand('seed',dsidx);
rp=randperm(numel(neginds));
rp=rp(1:min(numel(rp),5000));
negPats=ds.initPatches(neginds(rp),:);
initFeatsNeg=ds.initFeats(neginds(rp),:);

% loop over the detectors, initializing each one.  Note that ds.conf.params.graddescfun
% is called to actually do the initialization.
alldets=[];
featstokeep={};
for(i=1:size(ds.detectors{dsidx}.id,1))
  % Select a single detector.
  ctr=effstridx(ds.detectors{dsidx},i);
  disp(['optimizing: ' num2str(ctr.id)]);

  % Find the initial feature(s) for this patch (in this algorithm, it
  % will be the single initially sampled patch).
  data=[ds.initFeats(ds.initPatches(:,6)==myclusts(i),:);initFeatsNeg];

  % Compute whether each patch is a positive (it's associated with this detector)
  % or a negative (it came from initFeatsNeg).  
  lab=zeros(size(data,1),1);
  lastpos=sum(ds.initPatches(:,6)==myclusts(i));
  lab(1:lastpos)=1;

  % Normalize and rescale the detector.  TODO: this should probably
  % be a parameter.
  ctr.w=ctr.w/norm(ctr.w)*2;
  ctr.b=.1;

  % Call the gradient descent function to do the initialization.
  % doGradDescentProj will simply set ctr.b such that the constraint
  % is satisfied.
  [ctr_out_tmp,scores]=ds.conf.params.graddescfun(data',lab*2-1,[ctr.w ctr.b]',ones(size(lab')),0);
  ctr.w=c(ctr_out_tmp(1:end-1))';
  ctr.b=ctr_out_tmp(end);
  ctr_out(i)=ctr;

  % Read the scores produced by ds.conf.params.graddescfun for each patch
  % and select only the ones above -.02.  However, make sure that we have
  % at least numel(ctr.w)/5; too few patches may cause us to degenerate later.
  sscores=sort(scores,'descend');
  thr=min(-.02,sscores(min(ceil(numel(ctr.w)/5),numel(sscores))));
  dets=[posPats(ds.initPatches(:,6)==myclusts(i),:);negPats(scores((lastpos+1):end)>=thr,:)];
  featstokeep{end+1,1}=data([true(lastpos,1);(scores((lastpos+1):end)>=thr)'],:);
  dets(:,6)=myclusts(i);
  alldets=[alldets;dets];

end

% write out the initialized detectors & cache the detections in the 
% proper location.  If ds.sys.distproc.localdir is set, write the
% cached features there.  Note that this will "assign" this batch
% to one particular machine; any other machine that tries to access
% the cache will fail!
ds.round.detectors{dsidx}=str2effstr(ctr_out);
ds.round.prevdets{dsidx}=alldets;
if(dsfield(ds,'sys','distproc','localdir'))
  prevfeats=cell2mat(featstokeep);
  save([ds.sys.distproc.localdir 'prevfeats' num2str(dsidx) '_0.mat'],'prevfeats');
else
  ds.round.prevfeats{dsidx}=cell2mat(featstokeep);
end
