function res=det2mat(dets);
  if(isempty(dets))
    res=[];
    return;
  end
  dets=dets(:);
  dim=7;
  if(isfield(dets,'flip'))
    dim=8;
  end
  if(isfield(dets,'boxid'))
    dim=9;
  end

  res=zeros(numel(dets),dim);
  res(:,6)=[dets.detector];
  res(:,7)=[dets.imidx];
  if(isfield(dets,'flip'))
    res(:,8)=[dets.flip];
  end
  if(isfield(dets,'boxid'))
    res(:,9)=[dets.boxid];
  end
  pos=[dets.pos];
  res(:,1:5)=getBoxesForPedro(pos,[dets.decision]);
end
