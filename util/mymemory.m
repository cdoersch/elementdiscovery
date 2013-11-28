% Author: Carl Doersch (cdoersch at cs dot cmu dot edu)
%
% Print an estimate of memory usage for the current workspace.
%function res=mymemory
%ans=[];
memory_a=whos;
%if(nargout==1)
%else
  disp(['workspace memory usage: ' num2str(sum([memory_a.bytes])) ' bytes']);
  mem_total=num2str(sum([memory_a.bytes]));
%end
