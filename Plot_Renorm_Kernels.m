c = categorical({'EXP','PWL','QEXP','RAY'})

ax1 = subplot(4,3,1)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [0.2 0.1 0.6 0.0]
b = bar(c,score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [0 .5 .5]
b.EdgeColor = [0 .6 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['EXP  (', char(949), ' = 0.1)'])

ax2 = subplot(4,3,2)
c = categorical({'EXP','PWL','QEXP','RAY'});
score1 = [0.2 0.1 0.3 0.0]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 0 .5]
b.EdgeColor = [.6 0 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['EXP  (', char(949), ' = 0.01)'])

ax3 = subplot(4,3,3)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.2 0.1 0.3 0.0]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 .5 0]
b.EdgeColor = [.6 .6 0]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['EXP  (', char(949), ' = 0.001)'])

ax4 = subplot(4,3,4)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.2 0.3 1.0 0.3]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [0 .5 .5]
b.EdgeColor = [0 .6 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['PWL  (', char(949), ' = 0.1)'])

ax5 = subplot(4,3,5)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.2 0.3 1.0 0.3]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 0 .5]
b.EdgeColor = [.6 0 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['PWL  (', char(949), ' = 0.01)'])

ax6 = subplot(4,3,6)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.2 0.3 1.0 0.3]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 .5 0]
b.EdgeColor = [.6 .6 0]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['PWL  (', char(949), ' = 0.001)'])

ax7 = subplot(4,3,7)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.1 0.3 0.1 0.3]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [0 .5 .5]
b.EdgeColor = [0 .6 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['QEXP  (', char(949), ' = 0.1)'])

ax8 = subplot(4,3,8)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.1 0.3 0.1 0.3]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 0 .5]
b.EdgeColor = [.6 0 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['QEXP  (', char(949), ' = 0.01)'])

ax9 = subplot(4,3,9)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.1 0.4 0.1 0.3]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 .5 0]
b.EdgeColor = [.6 .6 0]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['QEXP  (', char(949), ' = 0.001)'])

ax10 = subplot(4,3,10)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.2 0.0 0.6 0.1]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [0 .5 .5]
b.EdgeColor = [0 .6 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['RAY  (', char(949), ' = 0.1)'])

ax11 = subplot(4,3,11)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.0 0.0 0.1 0.1]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 0 .5]
b.EdgeColor = [.6 0 .6]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['RAY  (', char(949), ' = 0.01)'])

ax12 = subplot(4,3,12)
c = categorical({'EXP','PWL','QEXP','RAY'})
score1 = [0.0 0.0 0.0 0.1]
b = bar(c, score1)%,'FaceColor',[0 .5 .5],'EdgeColor',[0 .6 .6],'LineWidth',1.5)
b.FaceColor = [.5 .5 0]
b.EdgeColor = [.6 .6 0]
b.LineWidth = 1.5
set(gca,'XTickLabel',{'EXP','PWL','QEXP','RAY'});
title(['RAY  (', char(949), ' = 0.001)'])
