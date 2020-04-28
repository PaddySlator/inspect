function [] = print_to_formats(figure_path,formats)
%envelope which prints the current figure (gcf) to different figure formats

%formats - cell of strings either 'fig' or handle for other image format e.g. '-depsc'

%remove whitespace from around figure
fig=gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


for i=1:length(formats)
    if strcmp(formats{i},'fig')
        savefig(figure_path);
    else
        print(gcf,figure_path,formats{i})
    end
end



end