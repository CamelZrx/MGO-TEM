function z = check_stop(eo,epsilon,delta)
% z is the exitflag of the algorithm, epsilon is the threshold for stopping
% of the experiment process, delta is the threshold for stopping of the
% sub-martingales.
global Su_MP MP CMP W1t Np mt UpP maxrun alpha_tt fmin
global Etrue_T plotflag num Up Etrue Ebest stop11 stop2a1

sfg = 0;                                   %stop flag for this check
if(length(CMP)>=1)    
    m = maxrun;                            %the number of equal value for sub-martingale to stop according to delta
    md = (1-exp(-1))/delta;
    stop1 = abs(MP(end)-CMP(end));         %stopping criteria 1 statistics
    if(length(W1t)<1)
        stop2a = 1;
    else
        stop2a = W1t(end);
    end
    
    stop11 = [stop11 stop1];
    stop2a1 = [stop2a1 stop2a];
    
    a = length(Su_MP);

    if(a>1)                                        %nonnegative sub-martingale has two or more elements
        alpha_t = Su_MP(end)-Su_MP(1)+exp(-1);     %up-crossing for the nonnegative sub-martingale    
        if(Su_MP(end)-Su_MP(end-1)==0)
            mt = [mt mt(end)+1];
        else
            mt = [mt 0];            
        end        
    elseif(a<=1)                            %nonnegative sub-martingale has one or less element
        alpha_t = exp(-1);
        mt = [mt 0];
    else                                       %both sub-martingales have two or more elements       
    end
    
    if(mt(end)>=m)
        sfg = 1;
    end
    
    UpP = 1;
    Etrue_T = [Etrue_T exp(-(eo-2*min(0,fmin)+1)/Np)];
    if(exp(-(eo-2*min(0,fmin)+1)/Np)>=alpha_t+eps)                  %upcrossing happens for nonnegative, check stopping criteria
        UpP = 0;
        if(stop1<=epsilon^2&&abs(stop2a)<delta&&mt(end)>=md)        %number, epsilon and sub-martingale converges
            sfg = 1;
        end
    end
        
    alpha_tt = [alpha_tt alpha_t+eps];
else
    alpha_tt = [alpha_tt exp(-1)+eps];
    UpP = 1;
    Etrue_T = [Etrue_T 0];
end

if(UpP==1)                  % do not upcross
    Up(num) = 1;
else
    Up(num) = 0;
end
    
z = sfg;

if(plotflag==1&&length(CMP)>=1)
    figure(1)
    clf
    delete(findall(gcf,'type','text'))
    
    sfh1=subplot(4,1,1);
    plot(stop11,'-kd','markersize',8,'MarkerFaceColor',[153 51 51]/255,...
        'MarkerEdgeColor',[0 0 0])
    hold on
    plot(epsilon^2*ones(1,length(MP(2:end)-CMP)),'-r')
    
    if num>=2        
        title(['Stopping Criteria 1 (\epsilon=' num2str(epsilon) ')'])
        xlim([1 max(2,length(MP(2:end)-CMP))])
        set(gca,'fontname','TimesNewRoman','fontsize',15)
        rx = xlim;
        ry = ylim;
        text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.6*(ry(2)-ry(1)),'Convergence Bound','fontname','TimesNewRoman','color','r','fontsize',15)
        x0=10;
        y0=205;
        width=900;
        height=400;
        set(gcf,'units','points','position',[x0,y0,width,height])
        sfh1.Position = sfh1.Position + [0 -0.02 0 0.03];
        axis auto
        xlim([1 max(2,length(MP(2:end)-CMP))])
        ybars = [epsilon^2 max(ylim)];
        patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1),...
        ybars(2) ybars(2)], [0.8 0.8 0.8],'FaceAlpha',.4,'Edgecolor','none')
    end
    
    sfh2 = subplot(4,1,2);
    bar(1:length(stop11),double(stop11<=epsilon^2),'EdgeColor','w','Barwidth',min(0.4*40/length(stop11),0.25))
    xlim([1 max(2,length(MP(2:end)-CMP))])
    if num>=2        
        title(['Stopping Criteria 1 Fired (1) or Not (0) (\epsilon=' num2str(epsilon) ')'])
        set(gca,'fontname','TimesNewRoman','fontsize',15)
        sfh2.Position = sfh2.Position + [0 0.05 0 -0.1];
        ylabel('S1')
        set(gca,'YTick',[0 1])
        rx = xlim;
    end
    ylim([-0.1 1.1])
    
    sfh3 = subplot(4,1,3);
    plot(abs(stop2a1),'-kd','markersize',8,'MarkerFaceColor',[153 51 51]/255,...
        'MarkerEdgeColor',[0 0 0])
    hold on
    plot(delta*ones(1,length(MP(2:end)-CMP)),'-r')
    if num>=2 
        rx = xlim;        
        title(['Stopping Criteria 2 (\delta=' num2str(delta) ')'])
        set(gca,'fontname','TimesNewRoman','fontsize',15)    
        sfh3.Position = sfh3.Position + [0 -0.02 0 0.03];
        text(rx(1)+0.1*(rx(2)-rx(1)),0.4,'Convergence Bound','fontname','TimesNewRoman','color','r','fontsize',15)
        ylim([-0.1 1])
    end
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ybars = [delta max(ylim)];
    patch([min(xlim) max(xlim) max(xlim) min(xlim)], [ybars(1) ybars(1),...
    ybars(2) ybars(2)], [0.8 0.8 0.8],'FaceAlpha',.4,'Edgecolor','none')
    
    sfh4=subplot(4,1,4);
    bar(1:length(stop2a1),double(abs(stop2a1)<=delta),'EdgeColor','w','Barwidth',min(0.4*40/length(stop2a1),0.25))
    ylim([-0.1 1.1])
    if num>=2        
        title(['Stopping Criteria 2 Fired (1) or Not (0) (\delta=' num2str(delta) ')'])
        set(gca,'fontname','TimesNewRoman','fontsize',15)
    end
    ylabel('S2')
    set(gca,'YTick',[0 1])
    sfh4.Position = sfh4.Position + [0 0.05 0 -0.1];
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ylim([-0.1 1.1])
    hold off
    drawnow
end

if(plotflag==1&&length(CMP)>=1)
    figure(2)
    clf    
    delete(findall(gcf,'type','text'))
    
    sfh1=subplot(4,1,1);
    plot(Etrue_T(2:end),'-kd','markersize',8,'MarkerFaceColor',[153 51 51]/255,...
        'MarkerEdgeColor',[0 0 0])
    hold on
    plot(alpha_tt(2:end),'-r')
    sfh1.Position = sfh1.Position + [0 -0.02 0 0.03];
    xlim([1 max(2,length(MP(2:end)-CMP))])
    rx = xlim;
    ry = ylim;
    text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.45*(ry(2)-ry(1)),'Transformed Local Minima'...
        ,'fontname','TimesNewRoman',...
        'color',[153 51 51]/255,'fontsize',15)
    text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.2*(ry(2)-ry(1)),'Convergence Bound'...
        ,'fontname','TimesNewRoman',...
        'color','r','fontsize',15)
    
    if num>=2        
        title('Stopping Criteria 3')
        set(gca,'fontname','TimesNewRoman','fontsize',15)
        x0=10;
        y0=205;
        width=900;
        height=400;
        set(gcf,'units','points','position',[x0,y0,width,height])
        axis auto
    end
    area(alpha_tt(2:end),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',.4,'Edgecolor','none')
    xlim([1 max(2,length(MP(2:end)-CMP))])
    
    sfh2=subplot(4,1,2);
    a = abs(Up-1);
    bar(1:length(a)-1,a(2:end),'EdgeColor','w','Barwidth',min(0.4*40/length(a),0.25))
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ylim([-0.1 1.1])
    if num>=2        
        title('Stopping Criteria 3 Fired (1) or Not (0)')
        set(gca,'fontname','TimesNewRoman','fontsize',15)
    end    
    ylabel('S3')
    set(gca,'YTick',[0 1])
    sfh2.Position = sfh2.Position + [0 0.05 0 -0.1];
    xlim([1 max(2,length(MP(2:end)-CMP))])
    
    sfh3=subplot(4,1,3);
    plot(mt(2:end),'-kd','markersize',8,'MarkerFaceColor',[153 51 51]/255,...
        'MarkerEdgeColor',[0 0 0])
    hold on
    plot((ceil(md)-1)*ones(1,length(MP(2:end)-CMP)),'-r')
    sfh3.Position = sfh3.Position + [0 -0.02 0 0.03];
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ylim([0 ceil(maxrun*1.4)])
    if num>=2        
        title(['Stopping Criteria 4a (\epsilon=' num2str(epsilon) ' and \delta=' num2str(delta) ')'])
        set(gca,'fontname','TimesNewRoman','fontsize',15)
        rx = xlim;
        ry = ylim;
        text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.85*(ry(2)-ry(1)),'Convergence Bound','fontname','TimesNewRoman','color','r','fontsize',15)
    end
    area((ceil(md)-1)*ones(1,length(MP(2:end)-CMP)),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',.4,'Edgecolor','none')

    sfh4=subplot(4,1,4);
    a = double(mt>=ceil(md));
    bar(1:length(a)-1,a(2:end),'EdgeColor','w','Barwidth',min(0.4*40/length(a),0.25))
    sfh4.Position = sfh4.Position + [0 0.05 0 -0.1];
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ylim([-0.1 1.1])
    if num>=2        
    title(['Stopping Criteria 4a Fired (1) or Not (0) (\delta=' num2str(delta) ')'])
    set(gca,'fontname','TimesNewRoman','fontsize',15)
    ylim([-0.1,1.1])
    end
    ylabel('S4a')
    set(gca,'YTick',[0 1])
    hold off
    drawnow
end

if(plotflag==1&&length(CMP)>=1)
    figure(3)
    clf
    delete(findall(gcf,'type','text'))
    
    sfh1=subplot(4,1,1);
    sfh1.Position = sfh1.Position + [0 -0.02 0 0.03];   
    plot(mt(2:end),'-kd','markersize',8,'MarkerFaceColor',[153 51 51]/255,...
        'MarkerEdgeColor',[0 0 0])
    hold on
    plot((ceil(maxrun)-1)*ones(1,length(MP(2:end)-CMP)),'-r')
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ylim([0 ceil(maxrun*1.4)])
    if num>=2        
        title(['Stopping Criteria 4b (\epsilon=' num2str(epsilon) ' and \delta=' num2str(delta) ')'])
        x0=10;
        y0=205;
        width=900;
        height=400;
        set(gcf,'units','points','position',[x0,y0,width,height])        
        set(gca,'fontname','TimesNewRoman','fontsize',15)
        rx = xlim;
        ry = ylim;
        text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.85*(ry(2)-ry(1)),'Convergence Bound','fontname','TimesNewRoman','color','r','fontsize',15)
    end
    xlim([1 max(2,length(MP(2:end)-CMP))])
    area((ceil(maxrun)-1)*ones(1,length(MP(2:end)-CMP)),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',.4,'Edgecolor','none')
    
    sfh2=subplot(4,1,2);
    a = double(mt>=ceil(maxrun));
    bar(1:length(a)-1,a(2:end),'EdgeColor','w','Barwidth',min(0.4*40/length(a),0.25))
    sfh2.Position = sfh2.Position + [0 0.05 0 -0.1];
    xlim([1 max(2,length(MP(2:end)-CMP))])
    ylim([-0.1 1.1])
    if num>=2        
    title(['Stopping Criteria 4b Fired (1) or Not (0) (\delta=' num2str(delta) ')'])
    set(gca,'fontname','TimesNewRoman','fontsize',15)
    ylim([-0.1,1.1])
    end
    ylabel('S4b')
    set(gca,'YTick',[0 1])
    hold off
    
    sfh3=subplot(4,1,3);
    l1 = [inf inf stop11]<=epsilon^2;
    l2 = [inf inf abs(stop2a1)]<=delta;
    l3 = [0 abs(Up-1)]>0;
    l4a = [0 mt]>=ceil(md);
    l4b = [0 mt]>=ceil(maxrun);
    a = double(l1&l2&l3&l4a|l4b);
    bar(1:length(a)-2,a(3:end),'EdgeColor','w','Barwidth',min(0.4*40/length(a),0.25))    
    sfh3.Position = sfh3.Position + [0 0.05 0 -0.1];
    if num>=2        
        title('((((S1&S2)&S3)&S4a)|S4b) Fired (1) or Not (0)')
        x0=10;
        y0=205;
        width=900;
        height=400;
        set(gcf,'units','points','position',[x0,y0,width,height])
        set(gca,'YTick',[0 1])
        axis auto
        set(gca,'fontname','TimesNewRoman','fontsize',15)
    end
    ylim([-0.1 1.1])
    xlim([1 max(2,length(MP(2:end)-CMP))])
    
    sfh4=subplot(4,1,4);
    plot(Etrue([1 3:end]),'-kd','markersize',8,'MarkerFaceColor',[153 51 51]/255,...
        'MarkerEdgeColor',[0 0 0])
    hold on
    plot(Ebest([1 3:end]),'-kd','markersize',8,'MarkerFaceColor',[0 51 102]/255,...
        'MarkerEdgeColor',[0 0 0])
    axis auto
    sfh4.Position = sfh4.Position + [0 0.01 0 0.03];
    rx = xlim;
    ry = ylim;
    text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.5*(ry(2)-ry(1)),'Found Local Minima'...
        ,'fontname','TimesNewRoman',...
        'color',[153 51 51]/255,'fontsize',15)
    text(rx(1)+0.1*(rx(2)-rx(1)),ry(1)+0.3*(ry(2)-ry(1)),'Running Minimum'...
        ,'fontname','TimesNewRoman',...
        'color',[0 51 102]/255,'fontsize',15)
    
    if num>=2        
        xlabel('Number of Samples')
        title('Optimization History')
        set(gca,'fontname','TimesNewRoman','fontsize',15)
        xlim([1 max(2,length(MP(2:end)-CMP))])
    end
end
num = num+1;
end