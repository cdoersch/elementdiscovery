function baseconf=overrideconf(baseconf,newconf)
  fnms=fieldnames(newconf)
  for(i=1:numel(fnms))
    if(isfield(baseconf,fnms{i})&&isstruct(getfield(baseconf,fnms{i}))&&isstruct(getfield(newconf,fnms{i})))
      baseconf=setfield(baseconf,fnms{i},overrideConf(getfield(baseconf,fnms{i}),getfield(newconf,fnms{i})));
    else
      baseconf=setfield(baseconf,fnms{i},getfield(newconf,fnms{i}));
    end
  end
end
