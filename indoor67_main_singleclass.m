% The same as indoor67_main.m, but it only runs the 'bowling'
% category.

global ds;
myaddpath;
finishedinit=false;
try,finishedinit=dsload('.ds.finishedinit');catch,end
if(~finishedinit)
  dscd('.ds');
  % Edit the following with your path to the indoor67 data.  This
  % directory should contain directories for airport_inside, artstudio,
  % etc.
  indoor67path='/code/data/indoor67/';
  % optionally, set the web-accessible path for the same directory--
  % i.e. the path such that [weburl '/airport_inside/{imagename}'] can
  % be used to download the image.  This is only used in displays;
  % if you just set it to the empty string, the code will still work,
  % but some links in the displays will be broken.
  weburl='';

  % set an output directory (on a filesystem that's accessible to all matlab's)
  dssetout('/PATH/TO/OUTPUT');
  % If you're running distributed and want to use SSH to send data between
  % nodes, use the following to set a directory where dswork can cache
  % files locally.  On Starcluster, /mnt/sgeadmin is a good place (assuming
  % you run your jobs as the user sgeadmin).
  %dssetlocaldir('/mnt/sgeadmin/');

  % reset the data store
  dsdelete('ds.*');

  % Set the number of parallel workers:
  njobs=12;
  % Set the amount of memory on each machine; dswork will estimate
  % how many jobs can be run simultaneously on each machine
  % based on this value.  Note that if you're running jobs alongside
  % the master, leave about 8GB for the master process
  ds.conf.mempermachine=59;
  % The machine where qsub can be run, or 'local' if you want to run on
  % a single machine:
  targmachine='local';

  % Any additional configuration to pass to dsmapredopen();
  distprocconf=struct();

  ds.masterscriptname=mfilename;

  if(~exist('dataset_indoor67_jpg.mat','file'))
    disp(['WARNING: you are about to create jpeg versions of the gifs in the indoor67 data directory.  '...
          'type ''return'' to continue or ''dbquit'' to cancel.']);
    keyboard
    preprocessindoor67(indoor67path,{'dataset_indoor67_jpg.mat','dataset_indoor67_test.mat'});
  end
  testdata=load('dataset_indoor67_test.mat');
  setdataset(20,testdata.imgs,indoor67path,testdata.labelnames,weburl);
  traindata=load('dataset_indoor67_jpg.mat');
  setdataset(19,traindata.imgs,indoor67path,traindata.labelnames,weburl);


  if(isfield(ds.conf.gbz{ds.conf.currimset},'imgsurl'))
    ds.imgsurl=ds.conf.gbz{ds.conf.currimset}.imgsurl;
  end
  s = RandStream('mcg16807','Seed',1234)
  RandStream.setDefaultStream(s) 
  ds.conf.params= struct( ...
    'ovlweight', 1, ... % use the inter-element communication scheme to set the weights.
    'negsperpos', 8, ... % during element training, the number of images we hard-mine 
                     ... % negatives from during for each positive training image.
    'maxpixels',300000,... % large images will be downsampled to this many pixels.
    'minpixels',300000,... % small images will be upsampled to this many pixels.
    'patchCanonicalSize', {[64 64]}, ... % resolution for each detected patch.
    'scaleIntervals', 8, ... % number of pyramid levels per scale.
    'sBins', 8, ... % pixel width/height for HOG cells
    'useColor', 1, ... % include Lab tiny images in the descriptor for a patch.
    'whitening', 1, ... % whiten the patch features
    'normbeforewhit', 1, ... % mean-subtract and normalize features before applying whitening
    'normalizefeats', 1, ... % mean-subtract and normalize features after applying whitening
    'graddescfun', @doGradDescentproj, ... % function ptr for the optimization function. It gets called
                                       ... % on each round of the optimization, including during
                                       ... % initialization.  See doGradDescentproj.m
    'stepsize', .05, ... % step size used by the optimizer
    'beta', 1, ... % beta value for the optimization
    'optimizationComputeLimit',1500,... % maximum number of vector-matrix multiplies that the
                                    ... % optimizer may perform on each training iteration
    'samplingOverlapThreshold', 0.6, ... % patches sampled initially can't have overlap larger
                                     ... % than this value.
    'samplingNumPerIm',20,... % sample this many patches per image.
    'multperim', 1, ... % allow multiple detections per image
    'nmsOverlapThreshold', 0.4 ... % overlap threshold for NMS during detection.
    )
  if(dsbool(ds.conf.params,'ovlweight'))
    ds.conf.params.lambdainit=.004;% lambda value for the optimization used during the first 
                        % training round (gets increased proportional to the number 
                        % of samples at later training rounds)
  else
    ds.conf.params.lambdainit=.02; % without the inter-element communication, everything
                                   % has a higher weight, and so we need to increase the
                                   % regularization.
  end

  classtorun='bowling';
  [~,classidxtorun]=ismember(classtorun,ds.conf.gbz{ds.conf.currimset}.labelnames);


  %pick which images to use out of the dataset
  imgs=ds.imgs{ds.conf.currimset};
  ds.myiminds=1:numel(imgs.label);

  % ds.roundinds stores the indices of the images to use at each iteration of training 
  [~,cls]=ismember(ds.imgs{ds.conf.currimset}.label(ds.myiminds),unique(ds.imgs{ds.conf.currimset}.label));
  imgsbyclass=distributeby(ds.myiminds(:),cls);
  rp=randperm(numel(ds.myiminds));
  % the first 3 rounds of training are just used to set the initial bandwidth, so we
  % use a very small subset of them.
  % actually it's pretty stupid that this takes 3 rounds; it could be done in 1 except
  % that doing so would generate tons of useless feature vectors in the current pipeline.
  ds.roundinds{1}=ds.myiminds(rp(1:20));
  ds.roundinds{2}=ds.myiminds(rp(21:60));
  ds.roundinds{3}=ds.myiminds(rp(61:120));
  % evenly divide the rest of the images. Note that the classes aren't quite balanced,
  % so the last round deals with extra/missing images.
  total_imgsbyclass = imgsbyclass;
  for(i=4:8)
    ds.roundinds{i}=[];
    for(j=1:numel(imgsbyclass))
      if(i==8)
        ul=numel(imgsbyclass{j});
      else
        ul=ceil(numel(total_imgsbyclass{j})/5);
      end
      ds.roundinds{i}=[ds.roundinds{i};imgsbyclass{j}(1:ul)];
      imgsbyclass{j}(1:ul)=[];
    end
    ds.roundinds{i}=ds.roundinds{i}(randperm(numel(ds.roundinds{i})));
  end
  dssave;
  if(~dsmapredisopen())
    dsmapredopen(njobs,targmachine,distprocconf);
    disp('waiting 10 sec for mapreducers to start...')
    pause(10)
  end

  % Generate the whitening matrix based on 1500 randomly sampled images.
  ds.aggcov.myiminds=ds.myiminds(rp(1:min(numel(rp),1500)));;
  dssave;
  dscd('ds.aggcov');
  dsrundistributed('aggregate_covariance',{'ds.myiminds'},struct('allatonce',1,'noloadresults',1,'maxperhost',max(1,floor(ds.conf.mempermachine/4))));
  total=0;
  clear featsum dotsum;
  dsload('ds.n');
  for(i=1:numel(ds.n))
    if(isempty(ds.n{i})),continue;end
    total=total+dsload(['ds.n{' num2str(i) '}'],'clear');
    if(~exist('dotsum','var'))
      dotsum=dsload(['ds.dotsum{' num2str(i) '}'],'clear');
    else
      dotsum=dotsum+dsload(['ds.dotsum{' num2str(i) '}'],'clear');
    end
    if(~exist('featsum','var'))
      featsum=dsload(['ds.featsum{' num2str(i) '}'],'clear');
    else
      featsum=featsum+dsload(['ds.featsum{' num2str(i) '}'],'clear');
    end
    if(any(isnan(dotsum(:)))||any(isnan(featsum(:))))
      keyboard;
    end
    disp(i);
  end
  covmat=(dotsum./total-(featsum'./total)*(featsum./total));
  covmat=covmat+.01*eye(size(covmat,1));
  dscd('.ds');
  ds.datamean=featsum./total;
  disp('performing matrix square root...');
  ds.invcovmat=inv(covmat);
  ds.whitenmat=sqrtm(ds.invcovmat);
  clear featsum dotsum total;
  clear covmat;
  dsdelete('ds.aggcov');
end
ds.finishedinit=true;
dssave;

initmodes=false;
try,initmodes=dsload('.ds.initmodes');catch,end
% This is the chunk of code that actually samples the patches that will become the 'modes',
% and sets up the positive/negative classes for each patch.
if(~initmodes)
  disp('sampling positive patches');
  dscd('.ds');
  setdataset(19);
  dsdelete('ds.round');
  dsdelete('ds.sample');
  dsdelete('ds.classperbatch');
  dsdelete('ds.detectors');
  dsdelete('ds.initFeats');
  dsdelete('ds.initPatches');
  dsdelete('ds.batchfordetr');
  dsdelete('ds.trainingrounds');
  dsdelete('ds.negprobperlabel');
  ulabels=unique(imgs.label);
  for(i=1:numel(ulabels))
    ds.negprobperlabel(ulabels(i))=ds.conf.params.negsperpos*sum(imgs.label(ds.myiminds)==ulabels(i))/sum(imgs.label(ds.myiminds)~=ulabels(i));
  end
  extrparams=ds.conf.params;
  ds.sample=struct();
  % Use this line to sample patches evenly from every image.
  % ds.sample.initInds=ds.myiminds;
  % alternatively, sample from just a few categories.  If you just want to sample candidate patches from a few categories,
  % note that you also need to sample some from negative images, just so that each patch detector will start with a few negatives.
  % I recommend around 4,000 negative patches to start off.
  positive_indices = ds.myiminds(imgs.label(ds.myiminds) == classidxtorun);
  negative_indices = ds.myiminds(imgs.label(ds.myiminds) ~= classidxtorun);
  rp=randperm(numel(negative_indices));
  ds.sample.initInds=[positive_indices, negative_indices(rp(1:200))];
  dsrundistributed('[ds.sample.patches{dsidx}, ds.sample.feats{dsidx}]=sampleRandomPatchesbb(ds.sample.initInds(dsidx),20);',{'ds.sample.initInds'},struct('maxperhost',max(1,floor(ds.conf.mempermachine/4))));

  % divide the sampled patches into batches (two images' worth of sampled
  % patches per batch).  The batches are purely for efficiency, specifically
  % to limit the number of files that get written.  Note that each
  % batch needs to have a single class label.
  batch_size=40;
  allpatches=cell2mat(ds.sample.patches(:));
  allfeats=cell2mat(ds.sample.feats(:));
  negfeats=[];
  negpatches=[];
  % optionally, get rid of any patches that you want to treat as purely
  % negative
  remove=ismember(allpatches(:,7),negative_indices);
  negpatches=allpatches(remove,:);
  allpatches(remove,:)=[];
  negfeats=allfeats(remove,:);
  allfeats(remove,:)=[];
  patchlabel=imgs.label(allpatches(:,7));
  ids=(1:size(allfeats,1))';
  [allpatches,allfeats,ids,patchlabel]=distributeby(allpatches,allfeats,ids,patchlabel);
  classperbatch={};
  for(i=1:numel(allpatches))
    batchpartition=c(bsxfun(@times,ones(batch_size,1),1:ceil(size(allpatches{i},1)/batch_size)));
    [allpatches{i},allfeats{i},ids{i}]=distributeby(allpatches{i},allfeats{i},ids{i},batchpartition(1:size(allpatches{i},1)));
    classperbatch{i}=repmat(patchlabel(i),size(allfeats{i},1),1);
  end
  allpatches=cat(1,allpatches{:});
  allfeats=cat(1,allfeats{:});
  ids=cat(1,ids{:});
  ds.classperbatch=cell2mat(classperbatch(:));
  initPatches=[cell2mat(allpatches);negpatches];
  disp(['sampled ' num2str(size(initPatches,1)) ' patches']);
  ds.initFeats=[cell2mat(allfeats);negfeats];

  % convert the patch features for each batch into a detector structure.
  ds.detectors=cellfun(@(x,y,z) struct('w',x,'b',y,'id',z),...
                     allfeats,...
                     cellfun(@(x) ones(size(x,1),1),allpatches,'UniformOutput',false),...
                     ids,...
                     'UniformOutput',false)';
  initPatches(1:numel(cell2mat(ids)),6)=(1:numel(cell2mat(ids)))';
  ds.initPatches=initPatches;
  % batchfordetr is an n-by-2 detector for the n detectors: column 1 is
  % a detector id, column 2 is the index of the batch containing it.
  ds.batchfordetr=[cell2mat(ids) cell2mat(cellfun(@(x,y) x*0+y,ids,c(num2cell(1:numel(ids))),'UniformOutput',false))];
  dssave();
  dsdelete('ds.sample')

  if(exist([ds.masterscriptname '_wait'],'file'))
    keyboard;
  end

  % initialize the set of detectors: this will only update the b value.
  ds.initFeats=[];
  runset=ds.sys.distproc.availslaves;
  dsrundistributed('autoclust_opt_init',{'ds.detectors'},struct('noloadresults',1,'maxperhost',max(1,floor(ds.conf.mempermachine/4)),'forcerunset',runset));
end
ds.initmodes=true;
dssave;

trainingrounds=false;
try,trainingrounds=dsload('.ds.trainingrounds');catch,end
if(~trainingrounds)
  dscd('.ds');
  setdataset(19);
  dsdelete('ds.gendpoolfeats');
  roundid=1;
  while(isfield(ds,['round' num2str(roundid)]))
    roundid=roundid+1;
  end
  uniquelabels=1:numel(ds.conf.gbz{ds.conf.currimset}.labelnames);
  ds.uniquelabels=uniquelabels(:)';
  while(roundid<=(numel(ds.roundinds)))
    thisroundmapreduce=false;
    try,thisroundmapreduce=dsload('.ds.round.thisroundmapreduce');catch,end
    if(~thisroundmapreduce)
      dsdelete('ds.round.myiminds');
      dsdelete('ds.round.roundid');
      dsdelete('ds.round.ndetrounds');
      dsdelete('ds.round.lambda');
      dsdelete('ds.round.beta');
      dsdelete('ds.round.newfeat');
      try
      dsdelete(['ds.round.newdets{' num2str(roundid) '}{1:' num2str(dssavestatesize('ds.round.newdets',2)) '}']);
      catch,end
      dsdelete('ds.nextround');
      ds.round.myiminds=ds.roundinds{roundid}; % images to run training on
      ds.round.ndetrounds=max(roundid-3,1); % the number of real detection rounds we've completed
      ds.round.roundid=roundid;
      if(~isfield(ds.round,'detrgroup'))
        % make a fake clustering for the first few rounds, where every patch gets its own cluster.
        ds.round.detrgroup=[ds.batchfordetr(:,1), (1:size(ds.batchfordetr,1))'];
      end
      if(roundid<=3)
        mph=max(1,floor(ds.conf.mempermachine/8));
      elseif(roundid<=4)
        mph=max(1,floor(ds.conf.mempermachine/2));
      else
        mph=max(1,floor(ds.conf.mempermachine/1.5));
      end
      if(mod(roundid,1)==0)
        % matlab's memory footprint grows even if it's not using the memory; restarting frees it.
        dsmapredrestart;
      end
      if(roundid>=4)
        % increase lambda (the bandwidth) proportional to the number of images we've run detection on.
        ds.round.lambda=(roundid-3)*ds.conf.params.lambdainit;
      else
        % if we're initializing, beta determines how large rho is.  Since we initialize on a small
        % set, we want to make beta artificially small so that the detector will fire less
        % when it starts doing detection for real.  This reduces the number of useless patch features
        % get generated.
        ds.round.beta=ds.conf.params.beta/3;
      end
      %end

      % the main mapreduce that runs detectors on images and sends the detected feature vectors
      % to reducers, each of which optimizes one batch of detectors. We use forcerunset to
      % make sure each detector is always optimized on the same machine, since these machines
      % cache features locally.
      dsmapreduce(['detectors=dsload(''ds.round.detectors'')'';'...
                   'imgs=dsload(''ds.imgs{ds.conf.currimset}'');'...
                   'dsload(''ds.classperbatch'');dsload(''ds.negprobperlabel'');'...
                   'posbats=find(imgs.label(ds.round.myiminds(dsidx))==ds.classperbatch);'... % run all detectors whose class matches the class of this image
                   'negbats=find(imgs.label(ds.round.myiminds(dsidx))~=ds.classperbatch);'...
                   'negbats=negbats(randsamplewithprob(ds.negprobperlabel(ds.classperbatch(negbats))));'... % run a random subset of the detectors for other classes
                   '[dets,feats]=detectInIm(effstrcell2mat(detectors([posbats(:); negbats(:)])),'...
                                            'ds.round.myiminds(dsidx),'...
                                            'struct(''thresh'',-.02/dsload(''ds.round.ndetrounds''),'...
                                            '''multperim'',dsload(''ds.round.roundid'')>2,'...
                                            '''flipall'',true));' ...
                   'ctridx=dsload(''ds.batchfordetr'');'...
                   'dsload(''ds.round.detrgroup'');'...
                   '[~,detrgroupord]=ismember(ds.round.detrgroup(:,1),ctridx(:,1));'...
                   'ovlweight=overlapReweightForImg(dets,[ctridx(:,1) ds.classperbatch(ctridx(:,2)) ds.round.detrgroup(detrgroupord,2)]);'...% ovlweights are the alpha_i,j from the paper
                   'ds.round.newfeat(1:numel(unique(ctridx(:,2))),dsidx)={struct(''assignedidx'',[],''feat'',[])};'...
                   'if(~isempty(dets)),'...
                     '[~,ctrpos]=ismember(dets(:,6),ctridx(:,1));'...
                     '[dets,feats,ovlweight,outpos]=distributeby(dets,single(feats),ovlweight,ctridx(ctrpos,2));'...
                     'ds.round.newfeat(outpos,dsidx)=cellfun(@(x,y,z) struct(''assignedidx'',x,''feat'',y,''ovlweights'',z),dets,feats,ovlweight,''UniformOutput'',false);'...
                   'end']...
                  ,'autoclust_optimize',{'ds.round.myiminds'},'ds.round.newfeat',struct('noloadresults',1,'forcerunset',runset),struct('maxperhost',mph),struct('maxperhost',max(1,floor(ds.conf.mempermachine/4))));
    end
    ds.round.thisroundmapreduce=true;
    dssave;
    % deleting prevfeats is optional (and prevents you from backing up the execution
    % to a previous round), but they take up a ton of space.
    dsdelete('ds.round.prevfeats');
    dsdelete('ds.round.component');
    dsdelete('ds.nextround.detrgroup');
    dsdelete('ds.nextround.reweights');
    dsdelete('ds.nextround.detsbyim');
    dsdelete('ds.nextround.prevweights');

    % Now that we've trained the detectors, we need to re-compute the
    % alpha values for each patch.
    if(roundid>=4)
      
      % findOverlapping3 finds overlapping clusters--i.e. it performs
      % the agglomerative element clustering step described in the paper.
      dsrundistributed(['dsload(''ds.classperbatch'');dsload(''ds.batchfordetr'');'...
                        'if(isempty(find(ds.classperbatch==ds.uniquelabels(dsidx)))),return;end,'...
                        '[~,~,ds.round.component{dsidx}]='...
                        'findOverlapping3(''ds.nextround.prevdets'',find(ds.classperbatch==ds.uniquelabels(dsidx)),'...
                        '[ds.batchfordetr(:,1),ds.classperbatch(ds.batchfordetr(:,2))],'...
                        'struct(''ndetsforoverlap'',.5,''maxoverlaps'',3,''clusterer'',''agglomerative''))'],'ds.uniquelabels',struct('maxperhost',max(1,floor(ds.conf.mempermachine/3))));

      % Construct the detrgroup matrix, which is an two-column matrix,  The left column
      % is the detector id, the right column specifies which cluster that element
      % belongs to.
      component=[];
      toadd=0;
      for(i=1:numel(ds.round.component))
        tmpcomponent=ds.round.component{i};
        if(numel(tmpcomponent)==0)
          continue;
        end
        tmpcomponent(:,2)=tmpcomponent(:,2)+toadd;
        component=[component;tmpcomponent];
        toadd=max(component(:,2));  
      end
      [~,cord]=ismember(ds.batchfordetr(:,1),component(:,1));
      component=component(cord,:);
      ds.nextround.detrgroup=component(:,1:2);

      % Collect the detections for each image and compute the alpha values.  Note
      % that if ds.conf.params.ovlweight is false or doesn't exist, then the
      % clustering generated above is ignored and weights are simply assigned
      % based on the number of detections for a given element in the image.
      detsbyim=cell2mat(dsload('ds.nextround.prevdets','clear')');
      [detsbyim,~,ord]=distributeby(detsbyim,detsbyim(:,7));
      ds.nextround.detsbyim=detsbyim';
      clear detsbyim;
      dsrundistributed(['ctridx=dsload(''ds.batchfordetr'');'...
                       'dsload(''ds.classperbatch'');'...
                       'dsload(''ds.nextround.detrgroup'');'...
                       '[~,detrgroupord]=ismember(ds.nextround.detrgroup(:,1),ctridx(:,1));'...
                       'ds.nextround.reweights{dsidx}=overlapReweightForImg(ds.nextround.detsbyim{dsidx},[ctridx(:,1) ds.classperbatch(ctridx(:,2)) ds.nextround.detrgroup(detrgroupord,2)]);'],'ds.nextround.detsbyim');
      ds.nextround.detsbyim={};
      dsload('ds.nextround.prevdets');
      reweights=mat2cell(invertdistributeby(ds.nextround.reweights(:),ord),cellfun(@(x) size(x,1),ds.nextround.prevdets),1);
      classfordetr(ds.batchfordetr(:,1))=ds.classperbatch(ds.batchfordetr(:,2));
      
      % The alpha value for the patch that initialized each cluster doesn't
      % follow the same rules as everything else, since we don't know when
      % the image it was sampled from will actually have detectors run
      % on it.  Hence, we instead start with a value of 1 and decay
      % exponentially toward the mean of the other weights.  I don't believe
      % this is important in the current implementation, but I haven't really
      % tried other approaches to setting this weight.
      for(i=1:numel(reweights))
        [currdets,currweights,detid,ord]=distributeby(ds.nextround.prevdets{i},reweights{i},ds.nextround.prevdets{i}(:,6));
        for(j=1:numel(currdets))
          ispos=find(ds.imgs{ds.conf.currimset}.label(currdets{j}(:,7))==c(classfordetr(currdets{j}(:,6))));
          if(ispos(1)~=1)
            error('detections got out of order');
          end
          if(numel(ispos)>1)
            currweights{j}(1)=.5^(roundid-4)+(1-.5^(roundid-4))*mean(currweights{j}(ispos(2:end)));
          end
        end
        reweights{i}=invertdistributeby(currweights,ord);
      end
      ds.nextround.prevweights=reweights(:)';
    end
    dsdelete(['ds.progressdisplay' num2str(roundid)]);

    % Generate a visualization of the progress.
    if(roundid>=4)
      % select a subset of the batches to display, since there's too many
      % detectors overall.
      batchestodisp=1:round(numel(unique(ds.batchfordetr(:,2)))/5):numel(unique(ds.batchfordetr(:,2)));
      batchestodisp=batchestodisp(1:min(10,numel(batchestodisp)));
      dets=cell2mat(dsload(['ds.nextround.prevdets{' num2str(batchestodisp) '}'])');
      ovlweights=cell2mat(dsload(['ds.nextround.prevweights{' num2str(batchestodisp) '}'])');
      [dets,ovlweights]=distributeby(dets,ovlweights,dets(:,6));
      for(i=1:numel(dets))
        [~,ord]=sort(dets{i}(:,5),'descend');
        ovlweights{i}=ovlweights{i}(ord(1:min(20,size(dets{i},1))));
        dets{i}=dets{i}(ord(1:min(20,size(dets{i},1))),:);
      end
      dets=cell2mat(dets);
      ovlweights=cell2mat(ovlweights);
      dsup(['ds.progressdisplay' num2str(roundid) '.patchimg'],extractpatches(double(dets)));
      conf=struct('dets',dets,'detrord',ds.batchfordetr(ismember(ds.batchfordetr(:,2),batchestodisp),1));
      if(dsbool(ds.conf.params,'ovlweight'))
        conf.ovlweights=ovlweights;
      end
      mhprender('patchdisplay.mhp',['ds.progressdisplay' num2str(roundid) '.displayhtml'],conf);
      fail=1;while(fail),try
      dssave;
      fail=0;catch ex,if(fail>5),rethrow(ex);end,fail=fail+1;end,end
      dsclear(['ds.progressdisplay' num2str(roundid)]);
    end
    ds.round=struct();
    dsmv('ds.round',['ds.round' num2str(roundid)]);
    %optionally, delete the rounds as we go along to save space.
    dsdelete(['ds.round' num2str(roundid) '.*']);
    dsmv('ds.nextround','ds.round');%keep the stub--it's how we measure progress!
    roundid=roundid+1;
  end
  ds.trainingrounds=true;
  dssave;
end

gendpoolfeats=false;
try,gendpoolfeats=dsload('.ds.gendpoolfeats');catch,end
if(~gendpoolfeats)
  dscd('.ds');
  dsdelete('ds.scores');
  dsdelete('ds.finids');
  dsdelete('ds.finmodel');
  dsdelete('ds.testpoolfeats');
  dsdelete('ds.poolfeats');
  dsdelete('ds.display_*');
  setdataset(19);
  % The element training is done.  Now we need to select the top elements
  % for each class.  Note that we perform the selection based on
  % ds.newdets, which we obtained by firing the detectors on images they
  % hadn't (yet) been trained on, so it's a roughly unbiased estimate
  % of their held out accuracy, but the scores from different rounds
  % of detections are not comparable.
  uniquelabels=1:numel(ds.conf.gbz{ds.conf.currimset}.labelnames);
  ds.uniquelabels=uniquelabels(:)';
  dsrundistributed(['dsload(''ds.batchfordetr'');dsload(''ds.classperbatch'');dsload(''ds.imgs'');dsload(''ds.roundinds'');'...
    'if(sum(ds.classperbatch==ds.uniquelabels(dsidx))==0),return;end,'...
    'dsload([''ds.newdets{'' num2str(numel(ds.roundinds)-2:numel(ds.roundinds)) ''}{'' num2str(find(ds.classperbatch(:)''==ds.uniquelabels(dsidx))) ''}'']);'...
    'detsbyround={};'...
    'for(i=numel(ds.roundinds)-1:numel(ds.roundinds)),'...
      'detsbyround{end+1,1}=structcell2mat(ds.newdets(i,:)'');'...
    'end,'...
    '[ds.finids{dsidx},ds.scores{dsidx}]=greedySelectDetrsCoverage(detsbyround,ds.imgs{ds.conf.currimset}.label==ds.uniquelabels(dsidx),.7,200,struct(''useoverlap'',1));'...
    'ds.newdets={}'...
    ],'ds.uniquelabels',struct('maxperhost',max(1,floor(ds.conf.mempermachine/4))));

  % Collect detections from the held out detectors and generate
  % a display for each class.  Note that this display isn't the same
  % as the one on the webpage.  This one displays held out detections
  % using the greedy ranking, whereas the one on the webpage displays
  % detections after training, only includes detections from the 
  % positive class, and ranks based on the SVM weight vectors that we haven't
  % even generated yet.
  heldoutdets={};
  for(i=numel(ds.roundinds)-2:numel(ds.roundinds))
    newdets=cell2mat(dsload(['ds.newdets{' num2str(i) '}{1:' num2str(dssavestatesize('ds.newdets',2)) '}'])');
    ds.newdets={};
    newdets(~ismember(newdets(:,6),cell2mat(ds.finids(:))),:)=[];
    heldoutdets{end+1,1}=newdets;
  end
  heldoutdets=cell2mat(heldoutdets(:));

  [heldoutdets,detid]=distributeby(heldoutdets,heldoutdets(:,6));
  for(i=1:numel(ds.finids))
    heldoutdetsbyclass=heldoutdets(ismember(detid,ds.finids{i}));
    ds.topk{i}=cell2mat(maxkall(heldoutdetsbyclass,5,20));
  end
  ds.classes=ds.conf.gbz{ds.conf.currimset}.labelnames(:)';
  dsrundistributed(['if(isempty(ds.finids{dsidx})),return;end,'...
                    'dsup([''ds.display_'' ds.classes{dsidx} ''.patchimg''],extractpatches(ds.topk{dsidx}));'...
                    'conf=struct(''dets'',ds.topk{dsidx},'...
                                '''detrord'',ds.finids{dsidx},'...
                                '''message'',{cellfun(@(x) [''score:'' num2str(x)],num2cell(ds.scores{dsidx}),''UniformOutput'',false)});'...
                    'mhprender(''patchdisplay.mhp'',[''ds.display_'' ds.classes{dsidx} ''.displayhtml''],conf);']...
                    ,{'ds.finids','ds.scores','ds.topk','ds.classes'},struct('noloadresults',true));

  % save out the final selected element detectors for easy access.
  model=effstrcell2mat(dsload('ds.round.detectors','clear')');
  model=selectdetrs2(model,cell2mat(ds.finids(:)));
  ds.finmodel=model;
  dssave;

  % generate a feature vector for each training image
  dsrundistributed('dsload(''ds.finmodel'');ds.poolfeats{dsidx}=distGenPooledFeats(ds.finmodel,ds.myiminds(dsidx))',{'ds.myiminds'},struct('noloadresults',true));
  trainlab=ds.imgs{ds.conf.currimset}.label;
  % switch to the testing set and
  setdataset(20);
  ds.mytestinds=1:numel(ds.imgs{ds.conf.currimset}.fullname);
  dsrundistributed('dsload(''ds.finmodel'');ds.testpoolfeats{dsidx}=distGenPooledFeats(ds.finmodel,ds.mytestinds(dsidx))',{'ds.mytestinds'},struct('noloadresults',true));
end
ds.gendpoolfeats=true;
dssave;
setdataset(20);
dscd('.ds');

%clear some memory
ds.round=struct();
dsmapredrestart;

% optionally load the ifv features generated by Junega, Vedaldi, Jawahar & 
% Zisserman 2013's implementation of ifv.  Note that that this takes about
% 5GB of extra memory, since we blindly load a 100,000-dimensional feature
% vector for every image.
if(0)
  % Note also we're assuming here that the ifv IMDB order is the same as
  % our order (since the ifv output doesn't actually include the IMDB :-/)
  % This will be true if the order of directory listings is deterministic,
  % since the IMDB gets their file listing in the same way we do.  However,
  % the madlab docs say that the order returned by dir actually depends on
  % the OS.
  ifvpermutation=1:10000;
  ifvpath='/PATH/TO/IFV';
  fils=cleandir([ifvpath '/data/codes/FKtest_comb_train_chunk*']);
  for(i=1:numel(fils))
    load([ifvpath '/data/codes/' fils(i).name]);
    [~,idx]=ismember(index,ifvpermutation);
    for(j=1:numel(idx))
      trainifvfeats(idx(j),:)=chunk(:,j)';
    end
    disp(i)
  end
  ifvkern=trainifvfeats*trainifvfeats';

  fils=cleandir([ifvpath '/data/codes/FKtest_comb_test_chunk*']);
  for(i=1:numel(fils))
    load([ifvpath '/data/codes/' fils(i).name]);
    [~,idx]=ismember(index,ifvpermutation);
    for(j=1:numel(idx))
      testifvFeats(idx(j),:)=chunk(:,j)';
    end
    disp(i)
  end
  ifvtestkern=trainifvfeats*testifvFeats';
  pfwt=.01;
  thresh=.6;
  ifvwt=18;
else
  ifvkern=0;
  ifvtestkern=0;
  ifvkern=0;
  ifvtestkern=0;
  if(dsbool(ds.conf.params,'ovlweight'))
    pfwt=.02;
    thresh=.5;
    ifvwt=0;
  else
    pfwt=.05;
    thresh=.3;
    ifvwt=0;
  end
end



trainpf2=cell2mat(dsload('ds.poolfeats','clear')');
testpf2=cell2mat(dsload('ds.testpoolfeats','clear')');
% the kernel function for the SVM.  It's actually a linear kernel after
% a rectified-linear feature transform.
kernfun=@(x,y) ((max(x+thresh,0))*(max(y+thresh,0))');
% the transformation that happens before the dot product in the kernel.
transfun=@(x) max(x+thresh,0);

% Everything is much more efficient if we work directly with the kernel matrix.
pfkern=kernfun(trainpf2,trainpf2);
kernmat=ifvkern*ifvwt+pfkern*pfwt;
classes=unique(trainlab)

% optionally you can open a matlabpool here, but if your'e using
% dswork your probably already have a lot of matlabs running on this machine.
% This honestly doesn't take that long, though.
% matlabpool open 12
parfor(i=1:numel(classes))
  inds=find(trainlab==classes(i));
  label=-ones(size(kernmat,1),1);
  label(inds)=1;
  trainedsvm{i}=svmtrain(label,double([(1:numel(label))' kernmat]),'-s 0 -t 4 -c .1 -h 0');
  disp(i)
end

% compute the scores for the testing images.
pftestkern=kernfun(trainpf2,testpf2);
traintestkern=ifvtestkern*ifvwt+pftestkern*pfwt;
testscr=[];
for(i=1:numel(trainedsvm))
  testscr(i,:)=(trainedsvm{i}.sv_coef'*traintestkern(trainedsvm{i}.SVs,:)-trainedsvm{i}.rho)*trainedsvm{i}.Label(1);
end

[~,label]=max(testscr,[],1);
truth=ds.imgs{ds.conf.currimset}.label';
% here we print the final performance.
perf=sum(truth==label)./numel(label)
dsdelete('ds.svm');
% Finally we generate heatmaps.
for(i=1:numel(classes))
  ds.svm{i}=getMinimalModel2(trainedsvm{i},transfun(trainpf2));
end
allscores=testscr';

dsdelete('ds.dispdets');
dsdelete('ds.easyimages');
dsdelete('ds.errorimages');
prevdets=cell2mat(dsload('ds.round.prevdets','clear')');
prevdets=distributeby(prevdets,prevdets(:,6));
ds.dispdets=cell2mat(maxkall(prevdets(:),5,5));
clear prevdets;
ds.dispdets=ds.dispdets(ismember(ds.dispdets(:,6),ds.finmodel.id),:);
dssave;
% there are two versions of the heatmap: first, we generate heatmaps for
% the most confident correctly-classified images, and then for the most confident
% errors. Note that the code below generates heatmaps for the detections
% that *lowered* the SVM score as well as those that raised it.  However, we
% omitted the ones that lowered the detection score from the website since
% they are harder to interpret.
for(disperror=[0 1])
  if(disperror)
    dscd('.ds.errorimages');
    errorval=allscores(sub2ind(size(allscores),(1:size(allscores,1))',label(:)))-allscores(sub2ind(size(allscores),(1:size(allscores,1))',truth(:)));
    ds.cls=[truth(:) label(:)]
  else
    dscd('.ds.easyimages');
    for(i=1:size(allscores,1))
      [scr,ds.cls(i,:)]=maxk(allscores(i,:),2);
      errorval(i)=scr(1)-scr(2);
    end
  end

  [confidence,ds.todisp]=maxk(errorval,100);
  ds.transfun=transfun;
  setdataset(20);

  dsrundistributed(['detrs=dsload(''.ds.finmodel'');svm=dsload(''.ds.svm'');transfun=dsload(''ds.transfun'');'...
                    'cls=dsload(''ds.cls'');dispdets=dsload(''.ds.dispdets'');dsload(''.ds.imgs'');'...
                    'im=im2double(getimg(ds.todisp(dsidx)));'...
                    'ds.origimg{dsidx}=im;'...
                    '[ds.leftposimg{dsidx},ds.leftnegimg{dsidx}]=dispClassifier('...
                      'detrs,im,svm{cls(ds.todisp(dsidx),1)},transfun,dispdets,[''ds.display_left_'' num2str(dsidx)],19);'...
                    '[ds.rightposimg{dsidx},ds.rightnegimg{dsidx}]=dispClassifier('...
                      'detrs,im,svm{cls(ds.todisp(dsidx),2)},transfun,dispdets,[''ds.display_right_'' num2str(dsidx)],19);'],'ds.todisp',struct('noloadresults',1));

  mhprender('errdisp.mhp','ds.errhtml',struct('trueclasses',{ds.conf.gbz{ds.conf.currimset}.labelnames(ds.cls(ds.todisp,1))},'guessclasses',{ds.conf.gbz{ds.conf.currimset}.labelnames(ds.cls(ds.todisp,2))},'confidence',confidence,'iserror',disperror));
  dssave;
end

dscd('.ds');
perf=sum(truth==label)./numel(label)
