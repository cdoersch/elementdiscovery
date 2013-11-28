ds_html{end+1}=sprintf('<table>\n');
ds_html{end+1}=sprintf('');
for(i=1:67)
ds_html{end+1}=sprintf('  <tr><td colspan="10"><h3 style="margin-top:20px;margin-bottom:2px;"><a stlye="dext-decoration:none;" href="../display_');
ds_html{end+1}=num2str([argv.classnames{i}]);
ds_html{end+1}=sprintf('/displayhtml.html">');
ds_html{end+1}=num2str([argv.classnames{i}]);
ds_html{end+1}=sprintf('</a><hr style="margin:0px;"/></td></tr>\n');
ds_html{end+1}=sprintf('  <tr>\n');
ds_html{end+1}=sprintf('  ');
for(j=1:10)
ds_html{end+1}=sprintf('<td><img src="patchimg[]/');
ds_html{end+1}=num2str([20*(i-1)+j]);
ds_html{end+1}=sprintf('.jpg"/></td>');
end
ds_html{end+1}=sprintf('  </tr>\n');
ds_html{end+1}=sprintf('  <tr>\n');
ds_html{end+1}=sprintf('  ');
for(j=11:20)
ds_html{end+1}=sprintf('<td><img src="patchimg[]/');
ds_html{end+1}=num2str([20*(i-1)+j]);
ds_html{end+1}=sprintf('.jpg"/></td>');
end
ds_html{end+1}=sprintf('  </tr>\n');
ds_html{end+1}=sprintf('\n');
ds_html{end+1}=sprintf('');
end
ds_reshtml=cell2mat(ds_html);
ds_html=[];
