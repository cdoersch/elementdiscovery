% [Aout, ..., uniquebins, inverse]=distributeby(A, ..., bins) will
% distribute A into a cell array, where rows of A that have the same
% values in the rows of bins will be mapped to the same cell of A.
%
% Specifically, uniquebins is the output of unique(bins,'rows'), where the
% rows are sorted, first by the first column, then by the second column,
% and so forth. For any index i, Aout{i} will be the (order-preserving) 
% concatenation of all A(j,:) where all(bins(j,:)==uniquebins(i,:)).  If A is 
% a struct, it is assumed to be a set of parallel arrays (i.e. an effstr); each 
% A{i} will be a struct with the same fields as the input, but with each 
% field containing only the values matching a particular value in bins.
% Finally, inverse is a cell array the same size as Aout, where inverse{i}(j)
% specifies which row of A corresponds to Aout{i}(j,:).  Note that
% invertdistributeby(Aout,inverse)=A.
%
% Additional arguments after A will have the exact same operation performed
% on them as A.
%
% [...]=distributeby(...,'descend') will invert the ordering of the 
% uniquebins output, and all the other variables as well.
function [varargout]=distributeby(varargin)
  %order='ascend';
  order=1;
  if(ischar(varargin{end})&&strcmp('descend',varargin{end}))
    %order='descend';
    order=-1;
    varargin(end)=[];
  end
  if(ischar(varargin{end})&&strcmp('ascend',varargin{end}))
    varargin(end)=[];
  end
  bins=varargin{end};
  mat=varargin(1:end-1);
  binid=sortrows(unique(bins,'rows'),(1:size(bins,2))*order);
  [~,idx]=ismember(bins,binid,'rows');
  [~,ord]=sort(idx);
  counts=histc(idx,1:size(binid,1));
  for(i=1:numel(mat))
    if(isstruct(mat{i}))
      fns=fieldnames(mat{i});
      restmp=cell(numel(counts),numel(fns));
      for(j=1:numel(fns))
        field=mat{i}.(fns{j});
        field=field(ord,:);
        restmp(:,j)=mat2cell(field,counts,size(field,2));
      end
      restmp=cell2struct(restmp,fns,2);
      varargout{i}=num2cell(restmp);
      %restmp=mat2cell(restmp,ones(numel(counts),1),size(restmp,2));
      %fn=fieldnames(mat{i});
      %restmp=cellfun(@(x) cell2struct(
    else
      res=mat{i}(ord,:);
      varargout{i}=mat2cell(res,counts,size(res,2));
    end
  end
  varargout{end+1}=binid;
  if(nargout>numel(varargout))
    varargout{end+1}=mat2cell(ord(:),counts,1);
  end
end
