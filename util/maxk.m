function [val,pos]=maxk(x,k)
  [val,pos]=sort(x,'descend');
  val=val(1:min(numel(x),k));
  pos=pos(1:min(numel(x),k));
end
