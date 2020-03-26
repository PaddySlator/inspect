function plot_3D_spectrum(F,grid,subplotdim,subplotind,gap,marg_h,marg_w)





% plot projections
%figure;hold on;

i=1;
%subtightplot(1,3,i,0.05,0.1, 0.05)
%subplot(1,3,i)
subtightplot(subplotdim(1),subplotdim(2),subplotind,gap,marg_h,marg_w) %,0.05,0.1, 0.05)
proj = squeeze(sum(F,1));
contour(grid{3},grid{2},proj)
%xlabel('T1 (ms)')
ylabel('T2* (ms)')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%axis square

i=2;
%subtightplot(1,3,i,0.05,0.1, 0.05)
%subplot(1,3,i) %,0.05,0.1, 0.05)
subtightplot(subplotdim(1),subplotdim(2),subplotind+1,gap,marg_h,marg_w)
proj = squeeze(sum(F,2));
contour(grid{3},grid{1},proj)
%xlabel('T1 (ms)')
ylabel('D (mm^2 s^{-1})')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%axis square


i=3;
%subtightplot(1,3,i,0.05,0.1, 0.05)
%subplot(1,3,i) %0.05,0.1, 0.05)
subtightplot(subplotdim(1),subplotdim(2),subplotind+2,gap,marg_h,marg_w)
proj = squeeze(sum(F,3));
contour(grid{2},grid{1},proj)
%xlabel('T2* (ms)')
ylabel('D (mm^2 s^{-1})')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%axis square


%hard-coded for my screen!
%set(gcf,'position',[-102        1119        1448         426])
%set(gcf,'position',[-102        1295         817         250])



end