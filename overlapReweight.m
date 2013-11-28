% Re-weight the detections specified by pos such that detections
% overlapping with many other detections receive low weight, and
% detections overlapping with few others receive high weight.
% The output weight for each patch (res2) is computed by (1) breaking
% the image into scale-space voxels, (2) assigning 'ownership' of each
% voxel like in the E-M algorithm, proportional to the weight of every
% patch covering the voxel (in the comments below, we call the 
% un-normalized ownership value for a patch the patch's 'claim' to a voxel),
% and (3) summing the ownership values across
% each patch.
function resw=overlapReweight(pos,origw,gridsize,minsize,conf)
  try
  if(size(pos,1)==0)
    resw=[];
    return;
  end
  pos=(pos-.5)/gridsize+1;
  minsize=prod(minsize/gridsize);
  imsize=ceil([max(pos(:,4)),max(pos(:,3))]);
  % group the detections such that detections with the same label do not
  % compete with each other.  Hence, in practice ownership for each
  % voxel is distributed among different detector groups rather than
  % among different detectors.
  if(isfield(conf,'groups'))
    [pos,origw,ord]=distributeby(pos,origw,(1:numel(origw))',conf.groups);
    ord=cell2mat(ord);
  else
    [pos,origw]=distributeby(pos,origw,(1:numel(origw))');
    ord=(1:numel(origw))';
  end
  imsize
  dim=[imsize ceil(log2(min(imsize))*5)];
  % field computes the sum of the weights the patches claiming each voxel.
  field=zeros(dim);
  field2=field;
  linidx=1;
  addvec={};
  for(j=1:numel(pos))
    allinds={};
    for(i=1:size(pos{j},1))
      p=pos{j}(i,:);
      patsize=(p(4)-p(2)+1)*(p(3)-p(1)+1);
      lvl=log2(sqrt(patsize/minsize))+1;
      p=[p([2 1]) lvl p([4 3]) lvl+5];
      f=floor(p);
      b=p-[f(1) f(2) f(3) f(1) f(2) f(3)];
      toadd=getaddvec(b).*origw{j}(i);
      if(size(pos{j},1)>1)
        % field2 is a temporary set of voxels used to compute one particular group's
        % claim to each voxel.
        field2(f(1):f(1)+size(toadd,1)-1,f(2):f(2)+size(toadd,2)-1,f(3):f(3)+size(toadd,3)-1)=max(field2(f(1):f(1)+size(toadd,1)-1,f(2):f(2)+size(toadd,2)-1,f(3):f(3)+size(toadd,3)-1),toadd);
        allinds{i}=[f(1),f(1)+size(toadd,1)-1,f(2),f(2)+size(toadd,2)-1,f(3), f(3)+size(toadd,3)-1];%inds;
      else
        try
        field(f(1):f(1)+size(toadd,1)-1,f(2):f(2)+size(toadd,2)-1,f(3):f(3)+size(toadd,3)-1)=field(f(1):f(1)+size(toadd,1)-1,f(2):f(2)+size(toadd,2)-1,f(3):f(3)+size(toadd,3)-1)+toadd;
        catch ex,dsprinterr;end
        linidx=linidx+1;
      end
    end
    if(size(pos{j},1)>1)
      field=field+field2;
      for(i=1:size(pos{j},1))
        f=allinds{i};
        % For each patch, keep track of the group's claim within the area of that patch.
        addvec{linidx}=c(field2(f(1):f(2),f(3):f(4),f(5):f(6)))./max(origw{j}(i),eps);%(fieldtmp(allinds{i})./max(origw{j}(i),eps));
        linidx=linidx+1;
      end
      field2(:)=0;
    end
  end
  % compute the inverse of the claims, since we compute the weight of a given patch as the sum over voxels
  % of a patch's claim to that voxel divided by the total claim to that voxel.
  field=1./max(field,eps);%bsxfun(@max,field,reshape(1/(prod(minsize)*2.^(0:size(field,3)-1),1,1,[])));

  pos=cell2mat(pos);
  origw=cell2mat(origw);
  resw=zeros(size(origw));
  % for each patch, re-compute its claim to each patch and multiply by
  % the values in field.
  for(i=1:size(pos,1))
    p=pos(i,:);
    patsize=(p(4)-p(2)+1)*(p(3)-p(1)+1);
    lvl=log2(sqrt(patsize/minsize))+1;
    p=[p([2 1]) lvl p([4 3]) lvl+5];
    f=floor(p);
    b=p-[f(1) f(2) f(3) f(1) f(2) f(3)];
      toadd=getaddvec(b);
      sta=sum(toadd(:));
      if(numel(addvec)>=i && ~isempty(addvec{i}))
        toadd=reshape(addvec{i},size(toadd));
      end
      if(sta>0)
        toadd=toadd./sta;
      end
      try
      resw(i)=resw(i)+(reshape(field(f(1):f(1)+size(toadd,1)-1,f(2):f(2)+size(toadd,2)-1,f(3):f(3)+size(toadd,3)-1),1,[])*toadd(:))*origw(i);
      catch ex,dsprinterr;end
    %end
  end
  if(any(resw>1.000001)),try,error('weightincrease');catch ex, dsprinterr;end,end
  if(any(isnan(resw))),try,error('nan');catch ex, dsprinterr;end,end
  resw=(min(resw,1).*origw);%
  resw(ord)=resw;
  if(any(resw<0)),try,error('less than zero');catch ex, dsprinterr;end,end
  catch ex,dsprinterr;end
end
%get the list of scale-space voxels for a particular detection, after it's been
%re-scaled into voxel coordinates.
function toadd=getaddvec(b)
      toadd=zeros(ceil(b(3:4)));
      dim=ceil(b(4:6));
      permvec=[];
      for(i=1:3)
        interpv=ones(dim(i),1);
        %interpv(1)=1-b(i);
        %interpv(end)=b(i+3)-numel(interpv)+1;
        permvec=[i permvec];
        if(i>1)
          interpv=permute(interpv,permvec);
        end
        mult{i}=interpv;
      end
      toadd=bsxfun(@times,bsxfun(@times,mult{1},mult{3}),mult{2});
      %toadd(2:end-1,2:end-1)=1;
      %toadd(1,2:end-1)=1-b(1);
      %toadd(end,2:end-1)=b(3)-size(toadd,2)+1;
      %toadd(2:end-1,1)=1-b(2);
      %toadd(2:end-1,end)=b(4)-size(toadd,2)+1;
      %toadd(1,1)=toadd(1,2)*toadd(2,1);
      %toadd(1,end)=toadd(2,end)*toadd(1,end-1);
      %toadd(end,1)=toadd(end,2)*toadd(end-1,1);
      %toadd(end,end)=toadd(end-1,end)*toadd(end,end-1);
end
