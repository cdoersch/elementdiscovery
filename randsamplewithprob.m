% sample from the indices 1:numel(probs) such that the probability of selecting i
% each one is approximately probs(i) (actually slightly larger due to rounding).
% Makes sure that res is nonempty.  
function res=randsamplewithprob(probs)
  expectednum=sum(probs)
  res=[];
  remainingprobs=expectednum;
  gotprobs=0;
  while(gotprobs<expectednum)
    idx=find(cumsum(probs)>rand*remainingprobs);
    remainingprobs=remainingprobs-probs(idx(1));
    gotprobs=gotprobs+probs(idx(1));
    probs(idx(1))=0;
    res=[res idx(1)];
  end
end
