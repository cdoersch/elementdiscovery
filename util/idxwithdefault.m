function res=idxwithdefault(data,inds,default)
  iscols=size(data,1)==0;
  data=[default;data(:)];
  inds(inds<=0)=0;
  inds=inds+1;
  res=data(inds);
  if(iscols)
    res=data';
  end
end
