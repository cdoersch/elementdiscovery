% Author: Carl Doersch (cdoersch at cs dot cmu dot edu)
%
% if you have a cell array where some cells contain structs of the
% same type and some cells are empty, matlab will crash.  This funciton
% gets rid of the empties so that cell2mat can function.
function res=structcell2mat(in)
  if(~iscell(in))
    res=in;
    return;
  end
  remove=cellfun(@(x) isempty(x),in);
  in(remove)=[];
  res=cell2mat(in);
end
