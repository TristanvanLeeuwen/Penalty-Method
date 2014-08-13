function savefig(figs,flabel)

tmp   = strfind(flabel,'/');
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

%close(figs);
%fclose all;