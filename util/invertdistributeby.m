% invert a distributeby operation, i.e. if you compute
% [Aout,~,inverse]=distributeby(A,bins)
% You can compute B as some function on each value in each cell of Aout
% and then call invertedistributeby(B,inverse) to get an array
% that's parallel with A again.
function val=invertdistributeby(val,ord)
  val=cell2mat(val);
  val(cell2mat(ord),:)=val;
end
