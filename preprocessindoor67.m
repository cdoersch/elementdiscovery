% Author: Carl Doersch (cdoersch at cs dot cmu dot edu)
%
% Generate an imgs structure and save it to gbz.datasetname.
% It reads all directories contained in gbz.cutoutdir and
% collects all of the images.
function preprocessindoor67(cutoutdir,datasetfiles)
for(imset=1:2)
  imdirs=dir(cutoutdir);
  [~,inds]=sort({imdirs.name})
  imdirs=imdirs(inds);
  imgs=[];
  imdata=[];
  smdata=[];
  if(imset==1)
    mitimgs=textread('TrainImages.txt','%s','delimiter','\n');
  else
    mitimgs=textread('TestImages.txt','%s','delimiter','\n');
  end
  for(fn=1:numel(imdirs))
    if(strcmp(imdirs(fn).name,'.')||strcmp(imdirs(fn).name,'..'))
      continue;
    end
    imdirs(fn).name
    imgs1=cleandir([cutoutdir '/' imdirs(fn).name]);
    [~,inds]=sort({imgs1.name});
    imgs1=imgs1(inds);
    rand('seed',fn);
    s=randperm(numel(imgs1));
    imgs2={};
    for(m=1:numel(imgs1))
      if(~ismember(lower([imdirs(fn).name filesep imgs1(m).name]),lower(mitimgs)))
        continue;
      end
      [~,pos]=ismember(lower([imdirs(fn).name filesep imgs1(m).name]),lower(mitimgs));
      mitimgs(pos)=[];
      if(dshassuffix(imgs1(m).name,'_gif.jpg'))
        fnm=[cutoutdir '/' imdirs(fn).name filesep imgs1(m).name];
        [im,map]=imread(fnm);
        if(size(map,2)>0)
          for(i=1:size(map,2))
            channel=map(:,i);
            data{i}=channel(im+1);
          end
          im=cat(3,data{:});
        end
        clear data
        fnm(end-6:end)='jpg.jpg';
        imgs1(m).name(end-6:end)='jpg.jpg';
        imwrite(im,fnm,'Quality',85);
        imgs1(m)=dir(fnm);
      end

      imgs2{end+1}.fullname=[imdirs(fn).name filesep imgs1(m).name];
      imtmp=imread([cutoutdir imgs2{end}.fullname]);
      sz=size(imtmp);
      imgs2{end}.imsize=sz(1:2);
      imgs2{end}.label=imdirs(fn).name;
    end
    imgs=[imgs;cell2mat(imgs2')];
    %imdata=[imdata; imdata1];
    %smdata=[smdata;smdata1];
  end
  imgs=str2effstr(imgs);
  labelnames=sort(unique(imgs.label));
  [~,imgs.label]=ismember(imgs.label,labelnames);
  save(datasetname{i},'imgs','labelnames');
end
