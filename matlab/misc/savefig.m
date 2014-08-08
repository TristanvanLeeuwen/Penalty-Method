function savefig(figs,flabel,n,m)

tmp=strfind(flabel,'/');
label = flabel(tmp(end)+1:end);

l = 1;
for k = figs
    figure(k);
    set(get(gca,'Xlabel'),'fontsize',20);
    set(get(gca,'Ylabel'),'fontsize',20)
    set(gca,'fontsize',20);
    try 
        set(get(gca,'Children'),'linewidth',2);
    catch
       % 
    end
    print(gcf,'-depsc',[flabel '_' num2str(char(96+l))]);
    l = l+1;
end
fid = fopen([flabel 'fig.tex'],'w');
fprintf(fid,'%%-----------------------------------------------------\n');
fprintf(fid,'\\newcommand{\\%s}[3]{\\begin{figure}[#1]\n',label);
fprintf(fid,'\\centering\n');
fprintf(fid,['\\begin{tabular}{' char(repmat(99,1,m)) '}\n']);
for l = 1:n
    fprintf(fid,['\\includegraphics[scale=#2]{%s.eps}'],[flabel '_' num2str(char(96-m+l*m+1))]);
    for k = 2:m
        fprintf(fid,['&\\includegraphics[scale=#2]{%s.eps}'],[flabel '_' num2str(char(96-m+l*m+k))]);
    end
    fprintf(fid,'\\\\\n');
    fprintf(fid,'{\\small (%s)}',num2str(char(96-m+l*m+1)));
    for k = 2:m
        fprintf(fid,'&{\\small (%s)}',num2str(char(96-m+l*m+k)));
    end
    fprintf(fid,'\\\\\n');
end

fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\caption{#3}\n');
fprintf(fid,'\\label{fig:%s}\n',label);
fprintf(fid,'\\end{figure}}\n');

%close(figs);
fclose all;