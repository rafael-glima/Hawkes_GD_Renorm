c = categorical({'EXP','PWL','QEXP','RAY'})

ax1 = subplot(1,4,1)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [ -9345.30 -9184.65 -4605.21 -4606.31;-8808.21 -8758.80 -4404.62 -4407.94; -6224.90 -5924.63 -3014.98 -3017.27;-9011.67 -9011.67 -4479.72 -4479.72]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
% b(3).EdgeColor = [0 .6 .6]
% b(4).EdgeColor = [.6 0 .6]
% b(5).EdgeColor = [0 .6 .6]
% b(6).EdgeColor = [.6 0 .6]
% b(7).EdgeColor = [0 .6 .6]
% b(8).EdgeColor = [.6 0 .6]4362.84
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['EXP  (', char(949), ' = 0.1)'])
l = cell(1,4);
l{1}='MLE'; l{2}=['RF-MLE (', char(949), ' = 0.1)']; l{3}=['RF-MLE (', char(949), ' = 0.01)']; l{4}=['RF-MLE (', char(949), ' = 0.001)'];
legend(b,l)

ax2 = subplot(1,4,2)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-6427.87 -6049.23 -6034.03 -6031.82; -8732.92 -8324.04 -8294.48 -8291.94; -5162.15 -3550.52 -3605.26 -3612.18; -6198.44 -6031.68 -6022.75 -6022.20]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
% b(3).EdgeColor = [0 .6 .6]
% b(4).EdgeColor = [.6 0 .6]
% b(5).EdgeColor = [0 .6 .6]
% b(6).EdgeColor = [.6 0 .6]
% b(7).EdgeColor = [0 .6 .6]
% b(8).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['PWL  (', char(949), ' = 0.1)'])

%%%% REVER VALORES PRA EPS0.01
ax3 = subplot(1,4,3)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-4263.04 -3577.93 -4249.27 -4263.04;-2246.73 -1623.26 -2246.73 -1994.41;-4834.63 -3824.03 -4824.63 -4834.63;-1898.07 -1898.07 -1898.07 -1898.07]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
% b(3).EdgeColor = [0 .6 .6]
% b(4).EdgeColor = [.6 0 .6]
% b(5).EdgeColor = [0 .6 .6]
% b(6).EdgeColor = [.6 0 .6]
% b(7).EdgeColor = [0 .6 .6]
% b(8).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['QEXP  (', char(949), ' = 0.1)'])

ax4 = subplot(1,4,4)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [-4726.30 -4726.30 -4726.30 -4726.30;-4439.19 -4439.19 -4439.19 -4439.19;-2991.98 -2931.48 -2985.19 -2991.98;-4507.61 -4507.61 -4507.61 -4507.61]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b(1).FaceColor = [0 .5 .5]
b(2).FaceColor = [.5 0 .5]
% b(3).FaceColor = [0 .5 .5]
% b(4).FaceColor = [.5 0 .5]
% b(5).FaceColor = [0 .5 .5]
% b(6).FaceColor = [.5 0 .5]
% b(7).FaceColor = [0 .5 .5]
% b(8).FaceColor = [.5 0 .5]
b(1).EdgeColor = [0 .6 .6]
b(2).EdgeColor = [.6 0 .6]
%b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['RAY  (', char(949), ' = 0.1)'])

