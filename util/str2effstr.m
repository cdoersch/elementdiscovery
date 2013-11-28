function res=str2effstr(str)
  nams=fieldnames(str);
  res=struct();
  for(i=1:numel(nams))
    if(numel(str)>0&&ischar(getfield(str(1),nams{i})))
      restruct={str.(nams{i})}';
    else
      restruct=cat(1,str.(nams{i}));
      if(isstruct(restruct))
        restruct=str2effstr(str);
      end
    end
    res=setfield(res,nams{i},restruct);
  end
end
