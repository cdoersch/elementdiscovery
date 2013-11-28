function model=selectdetrs2(model,ovl)
  [~,inds]=ismember(ovl,[model.id]);
  %model.w=model.w(inds,:);
  %model.rho=model.rho(inds);
  %model.id=model.id(inds);
  model=effstridx(model,inds);
end
