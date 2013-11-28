function data=maxkall(data,idx,nmax)
try
  for(i=1:numel(data))
    [~,ord]=mink(data{i}(:,idx),nmax);
    data{i}=data{i}(ord,:);
  end
catch ex,dsprinterr;end
end
