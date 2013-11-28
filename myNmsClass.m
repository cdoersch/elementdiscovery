function pick=myNmsClass(dets,overlap)
  [dets,~,ord]=distributeby(dets(:,1:5),dets(:,6));
  for(i=1:numel(dets))
    pick{i,1}=ord{i}(myNms(dets{i},overlap));
  end
  pick=cell2mat(pick);
end
