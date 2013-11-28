function res=genheatmap(wt,pos,sz)
  res=zeros(sz(1),sz(2));
  pos=round(pos);
  for(i=1:numel(wt))
    res(pos(i,2):pos(i,4),pos(i,1):pos(i,3))=res(pos(i,2):pos(i,4),pos(i,1):pos(i,3))+wt(i)/((pos(i,3)-pos(i,1)+1)*(pos(i,4)-pos(i,2)+1));
  end
end
