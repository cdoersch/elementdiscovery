% Author: Carl Doersch (cdoersch at cs dot cmu dot edu)
%
% Extract patches corresponding to detections in images.
% detsimple is the set of detections, with a 'pos' field
% specifying positions, and a 'imidx' field, specifying
% an index into the current image set. imgs_all us unused.
function res=extractpatches(detsimple,imgs_all,conf)
global ds;
if(~isstruct(detsimple))
  detsimple=mat2det(detsimple);
end
if(~exist('conf','var'))
  conf=struct();
end
try
%    dsload('ds.imgs');
%    global imgs;
%    if(isempty(imgs))
  tmpimgs=dsload('.ds.imgs{ds.conf.currimset}');
  imgs=cell(numel(tmpimgs.fullname),1);
%    end
numel(detsimple)
loaded=[];
ct=1;
if(~dsbool(conf,'featurize'))
  res=cell(1,numel(detsimple));
end
[~,ord]=sort([detsimple.imidx]);
  for(dsidx=ord(:)')
    pos=detsimple(dsidx).pos;
    i=detsimple(dsidx).imidx;
    if(i==0)
      if(isfield(conf,'explicitpatch')&&~isempty(conf.explicitpatch{dsidx}))
        res{dsidx}=conf.explicitpatch{dsidx};
      else
        res{dsidx}=cat(3,1,1,1);
      end
      if(~dsbool(conf,'noresize'))
        pat=imresize(res{dsidx},[reszy reszx]);
      end
      continue;
    end
    if(numel(imgs)<i||isempty(imgs{i}))
      imgs{i}=getimg(ds,i);%imread([ds.conf.gbz.cutoutdir imgs_all(i).fullname]);
      currinds=find([detsimple(ord).imidx]==i);
      postmp=[detsimple(ord(currinds)).pos];
      padx=max(0,max(-min([postmp.x1])+1,max([postmp.x2]-size(imgs{i},2))));
      pady=max(0,max(-min([postmp.y1])+1,max([postmp.y2]-size(imgs{i},1))));
      if(padx>0||pady>0)
        imgs{i}=padarray(imgs{i},[pady padx 0],'replicate');
      end
      loaded(loaded==i)=[];
      loaded=[loaded i];
      if(numel(loaded)>1)
        imgs{loaded(1)}=[];
        loaded(1)=[];
      end
    end
    %if(padsize>0)
      %toaddy=round((pos.y2-pos.y1+1)*padsize);
      %toaddx=round((pos.x2-pos.x1+1)*padsize);
      pos.y1=pos.y1+pady;
      pos.y2=pos.y2+pady;
      pos.x1=pos.x1+padx;
      pos.x2=pos.x2+padx;
    %end
    if(dsbool(conf,'noresize'))
      pat=imgs{i}(pos.y1:pos.y2,pos.x1:pos.x2,:);
    else
      maxsz=max(pos.y2-pos.y1,pos.x2-pos.x1);
      if(dsbool(conf,'explicitsize'))
        reszx=conf.explicitsize(2);
        reszy=conf.explicitsize(1);
      else
        reszx=80;%*(pos.x2-pos.x1)/maxsz;
        reszy=80;%*(pos.y2-pos.y1)/maxsz;
      end
      pat=imresize(imgs{i}(pos.y1:pos.y2,pos.x1:pos.x2,:),[reszy reszx]);
    end
    if(dsbool(detsimple(dsidx),'flip'))
      pat=pat(:,end:-1:1,:);
    end
    if(dsbool(conf,'featurize'))
      conf.imid=i;
      feat=patch2feat({pat},conf);
      if(~exist('res','var'))
        res=zeros(numel(detsimple),numel(feat),'single');
        res(dsidx,:)=single(feat);
      else
        res(dsidx,:)=single(feat);
      end
    else
      res{dsidx}=pat;
    end
    if(mod(ct,10)==0)
      disp(ct)
    end
    ct=ct+1;
    
  end
  dsclear('.ds.overalldets');
catch ex
  dsprinterr
end
