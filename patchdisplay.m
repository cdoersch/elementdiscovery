ds_html{end+1}=sprintf('<html><body>\n');
ds_html{end+1}=sprintf('');
if(isfield(argv,'mainimgs'))
ds_html{end+1}=sprintf('  <table><tr>\n');
ds_html{end+1}=sprintf('    ');
for(i=1:numel(argv.mainimgs))
ds_html{end+1}=sprintf('      <td><img src="');
ds_html{end+1}=num2str([argv.mainimgs{i}]);
ds_html{end+1}=sprintf('"/></td>\n');
ds_html{end+1}=sprintf('    ');
end
ds_html{end+1}=sprintf('  </tr></table>\n');
ds_html{end+1}=sprintf('');
end
ds_html{end+1}=sprintf('<table>\n');
ds_html{end+1}=sprintf('');
 
if(~isfield(argv,'ovlweights'))
  argv.ovlweights=zeros(size(argv.dets(:,[])));
end
[dets,posinpatchimg,ovlweights,detid]=distributeby(argv.dets,(1:size(argv.dets,1))',argv.ovlweights,argv.dets(:,6));
if(isfield(argv,'detrord'))
  [~,idxord]=ismember(argv.detrord,detid);
else
  idxord=1:numel(detid);
  argv.detrord=detid;
end
if(~isfield(argv,'message'))
  argv.message=repmat({''},numel(argv.detrord),1);
end
gbz=dsload('ds.conf.gbz{ds.conf.currimset}');
imgs=dsload('.ds.imgs{ds.conf.currimset}');
if(~isfield(gbz,'imgsurl'))
  gbz.imgsurl='';
end
for(i=1:numel(idxord)) 
ds_html{end+1}=sprintf('  <tr><td>detector no. ');
ds_html{end+1}=num2str([argv.detrord(i)]);
ds_html{end+1}=sprintf(' <br/> ');
ds_html{end+1}=num2str([argv.message{i}]);
ds_html{end+1}=sprintf('</td>\n');
ds_html{end+1}=sprintf('  ');
if(idxord(i) ~= 0)
    curdets=dets{idxord(i)};
    curpos=posinpatchimg{idxord(i)};
    curwt=ovlweights{idxord(i)};
    [~,ord]=sort(curdets(:,5),'descend');
    for(j=1:size(curdets,1)) 
ds_html{end+1}=sprintf('      <td>\n');
ds_html{end+1}=sprintf('        <a href="');
if(isfield(argv,'url'))
ds_html{end+1}=sprintf('');
ds_html{end+1}=num2str([argv.url{curpos(ord(j))}]);
ds_html{end+1}=sprintf('');
else
ds_html{end+1}=sprintf('');
ds_html{end+1}=num2str([gbz.imgsurl]);
ds_html{end+1}=sprintf('/');
ds_html{end+1}=num2str([imgs.fullname{curdets(ord(j),7)}]);
ds_html{end+1}=sprintf('');
end
ds_html{end+1}=sprintf('">\n');
ds_html{end+1}=sprintf('          <img src="patchimg[]/');
ds_html{end+1}=num2str([curpos(ord(j))]);
ds_html{end+1}=sprintf('.jpg" title="');
ds_html{end+1}=num2str([curdets(ord(j),5)]);
ds_html{end+1}=sprintf(' ');
ds_html{end+1}=num2str([curwt(ord(j),:)]);
ds_html{end+1}=sprintf('"/>\n');
ds_html{end+1}=sprintf('        </a>\n');
ds_html{end+1}=sprintf('      </td>\n');
ds_html{end+1}=sprintf('    ');
end
  end
ds_html{end+1}=sprintf('  </tr>\n');
ds_html{end+1}=sprintf('');
end
ds_html{end+1}=sprintf('      \n');
ds_html{end+1}=sprintf('</table>\n');
ds_html{end+1}=sprintf('</body></html>\n');
ds_reshtml=cell2mat(ds_html);
ds_html=[];
