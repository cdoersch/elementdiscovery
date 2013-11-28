% set the dataset that the algorithm uses.  In particular, this
% affects the behavior of getimgs and getimg.  Note that this sets
% the image set only for the current directory that dswork is dscd'd
% into.  
%
% idx is the index of the image set in the .ds.imgs array.  It should
% be unique for every dataset you use.  
%
% imgs is a struct with the fields:
%   - 'fullname': a cell array of labels, such that 
%     [datadir '/' imgs.fullname{i}] will resolve to the  absolute 
%     path of an image.
%   - 'label': an index in labelnames for the image's label.
%
% datadir: the directory containing all of your data
%
% labelnames: names for each label used in the dataset.
%
% weburl: an optional url that points to the web-accessible
% location of datadir, i.e. such that [imgsurl '/' imgs.fullname{i}]
% lets you download the image.  Used in html displays.
function setdataset(idx,imgs,datadir,labelnames, weburl)
  global ds;
  dsup('ds.conf.currimset',idx);
  if(nargin>1)
    if(~exist('weburl','var'))
      weburl='';
    end
    gbz=struct('labelnames',{labelnames},'imgsurl',weburl,'cutoutdir',datadir);
    mydir=dspwd;
    dscd('.ds');
    dsup(['ds.imgs{' num2str(idx) '}'],imgs);
    dsup(['ds.conf.gbz{' num2str(idx) '}'],gbz);
    dscd(mydir);
    dsup(['ds.conf.gbz{' num2str(idx) '}'],gbz);% TODO: this should really only be kept in 
                                                % the root .ds.conf; relying on the
                                                % non-root ds.conf creates
                                                % unintuitive behavior.
  end
end
