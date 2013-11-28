function res=mat2det(dets);
  res.decision=dets(:,5);
  res.pos.x1=dets(:,1);
  res.pos.x2=dets(:,3);
  res.pos.y1=dets(:,2);
  res.pos.y2=dets(:,4);
  res.imidx=dets(:,7);
  res.detector=dets(:,6);
  if(size(dets,2)>=8)
    res.flip=dets(:,8);
  end
  if(size(dets,2)>=9)
    res.boxid=dets(:,9);
  end
  res=effstr2str(res);
end
