% efficient nearest-neighbors in Euclidean distance.
% each row of toassign is assigned to the nearest row in targets.
% closest(i) is the row-index in targets of the closest element 
% for toassign(i,:).  outdist(i) is the distance to that point.
function [closest,outdist]=assigntoclosest(toassign,targets,nonrm)
  global ds;
  if(isempty(targets))
    closest=[];
    outdist=[];
    return;
  end
  targsq=targets.^2;%sum(targets.^2,2);
  closest=zeros(size(toassign,1),1);
  outdist=zeros(size(toassign,1),1);
  for(i=1:800:size(toassign,1))
    inds=i:min(i+800-1,size(toassign,1));
    batch=toassign(inds,:);
    batchsq=sum(batch.^2,2);
    inprod=targets*(batch');
    if(dsbool(ds.conf,'whiteningv2')||(exist('nonrm','var')&&nonrm))
      dist=inprod;
      [outdist(inds),closest(inds)]=max(dist,[],1);
    else
      %dist=bsxfun(@plus,bsxfun(@minus,batchsq',2*inprod),targsq);
      normval=sqrt(bsxfun(@rdivide,targsq*(batch'~=0),sum(batch'~=0,1))-bsxfun(@rdivide,(targets*(batch'~=0)).^2,sum(batch'~=0,1).^2));
      %normval=sqrt(targsq*(batch'~=0)-bsxfun(@rdivide,(targets*(batch'~=0)).^2,sum(batch'~=0,1)));

      normval(normval==0)=1;
      dist=(-bsxfun(@rdivide,inprod,sum(batch'~=0,1))./normval);
      %if(any(dist(:))<0)
        %keyboard;
      %end
      [outdist(inds),closest(inds)]=min(dist,[],1);
    end
  end
end
