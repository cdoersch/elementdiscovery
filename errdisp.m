ds_html{end+1}=sprintf('<html><body>\n');
ds_html{end+1}=sprintf('<table>\n');
ds_html{end+1}=sprintf('  ');
for(i=1:numel(ds.todisp)) 
ds_html{end+1}=sprintf('    <tr>\n');
ds_html{end+1}=sprintf('      <td><img src="');
ds_idxstr=i;
ds_html{end+1}=dsreldiskpath(['ds.origimg' '{' num2str(ds_idxstr) '}'],outdirpath);
if(numel(ds_html{end})>0),ds_html{end}=ds_html{end}{1};else,ds_html{end}='';end
ds_html{end+1}=sprintf('"/></td>\n');
ds_html{end+1}=sprintf('      <td><a href="');
ds_html{end+1}=num2str(['display_left_' num2str(i) '_pos/displayhtml.html']);
ds_html{end+1}=sprintf('"><img src="');
ds_idxstr=i;
ds_html{end+1}=dsreldiskpath(['ds.leftposimg' '{' num2str(ds_idxstr) '}'],outdirpath);
if(numel(ds_html{end})>0),ds_html{end}=ds_html{end}{1};else,ds_html{end}='';end
ds_html{end+1}=sprintf('"/></a><br/>');
if(~argv.iserror)
ds_html{end+1}=sprintf('top guess');
else
ds_html{end+1}=sprintf('ground truth');
end
ds_html{end+1}=sprintf(' (pos) (confindence ');
ds_html{end+1}=num2str([argv.confidence(i)]);
ds_html{end+1}=sprintf('):');
ds_html{end+1}=num2str([argv.trueclasses{i}]);
ds_html{end+1}=sprintf('</td>\n');
ds_html{end+1}=sprintf('      <td><a href="');
ds_html{end+1}=num2str(['display_left_' num2str(i) '_neg/displayhtml.html']);
ds_html{end+1}=sprintf('"><img src="');
ds_idxstr=i;
ds_html{end+1}=dsreldiskpath(['ds.leftnegimg' '{' num2str(ds_idxstr) '}'],outdirpath);
if(numel(ds_html{end})>0),ds_html{end}=ds_html{end}{1};else,ds_html{end}='';end
ds_html{end+1}=sprintf('"/></a><br/>');
if(~argv.iserror)
ds_html{end+1}=sprintf('top guess');
else
ds_html{end+1}=sprintf('ground truth');
end
ds_html{end+1}=sprintf(' (neg):');
ds_html{end+1}=num2str([argv.trueclasses{i}]);
ds_html{end+1}=sprintf('</td>\n');
ds_html{end+1}=sprintf('      <td><a href="');
ds_html{end+1}=num2str(['display_right_' num2str(i) '_pos/displayhtml.html']);
ds_html{end+1}=sprintf('"><img src="');
ds_idxstr=i;
ds_html{end+1}=dsreldiskpath(['ds.rightposimg' '{' num2str(ds_idxstr) '}'],outdirpath);
if(numel(ds_html{end})>0),ds_html{end}=ds_html{end}{1};else,ds_html{end}='';end
ds_html{end+1}=sprintf('"/></a><br/>');
if(~argv.iserror)
ds_html{end+1}=sprintf('second guess');
else
ds_html{end+1}=sprintf('guess');
end
ds_html{end+1}=sprintf(' (pos):');
ds_html{end+1}=num2str([argv.guessclasses{i}]);
ds_html{end+1}=sprintf('</td>\n');
ds_html{end+1}=sprintf('      <td><a href="');
ds_html{end+1}=num2str(['display_right_' num2str(i) '_neg/displayhtml.html']);
ds_html{end+1}=sprintf('"><img src="');
ds_idxstr=i;
ds_html{end+1}=dsreldiskpath(['ds.rightnegimg' '{' num2str(ds_idxstr) '}'],outdirpath);
if(numel(ds_html{end})>0),ds_html{end}=ds_html{end}{1};else,ds_html{end}='';end
ds_html{end+1}=sprintf('"/></a><br/>');
if(~argv.iserror)
ds_html{end+1}=sprintf('second guess');
else
ds_html{end+1}=sprintf('guess');
end
ds_html{end+1}=sprintf(' (neg):');
ds_html{end+1}=num2str([argv.guessclasses{i}]);
ds_html{end+1}=sprintf('</td>\n');
ds_html{end+1}=sprintf('    </tr>\n');
ds_html{end+1}=sprintf('  ');
end 
ds_html{end+1}=sprintf('</table>\n');
ds_html{end+1}=sprintf('</body></html>\n');
ds_reshtml=cell2mat(ds_html);
ds_html=[];
