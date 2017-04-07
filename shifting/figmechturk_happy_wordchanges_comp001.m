function figmechturk_happy_wordchanges_comp001(bigindex,freq1,freq2,file_tag, ...
                                        title_text,num);
% figmechturk_happy_wordchanges_comp001(bigindex,freq1,freq2,file_tag,title_text,num);
%
% works with mechanical turk word list
%
% prints a valence shift word graph for freq2 with respect to freq1
% placing most positive text on right
% 
% outfile is named figmechturk_happy_wordchanges_comp001/figmechturk_happy_wordchanges_comp001_[file_tag]_noname.pdf
%
% num is the number of shifts to show in the bar graph
%
% title_text should be length 2 cell containing strings describing
% each distribution in brief

picdir = '/Users/danforth/Pictures/science/twitter/figures/danforth/paper-figs';

more off;

figure(1);
clf;
figwidth = 800;
figheight = 1350;
figshape(figwidth,figheight);
aspectratio = figheight/figwidth;

mainaxeswidth = 0.5;
mainaxeslength = 0.76*num/50;

insertwidth = 0.11;
boxplotlength = 0.04;
circleplotlength = 0.08;
cumulplotlength = 0.15;

%% boxplotlength = insertwidth/2;
%% circleplotlength = insertwidth;

%% verticalseplength = 0.025;
verticalseplength = 0.02;
sideaxeslength =  mainaxeslength - boxplotlength - circleplotlength ...
    - 2*verticalseplength;

%% main plot
%% vertical axis length
axesinfo(1).position = [.1, .12, mainaxeswidth, mainaxeslength];
%axesinfo(1).position = [.1, .12, mainaxeswidth, mainaxeslength];

%% relative text size (number of mech turk words)
axesinfo(2).position = [.49, .12 + circleplotlength + 3.5*verticalseplength, insertwidth, boxplotlength];

%% four circles
axesinfo(3).position = [.485, .12 + verticalseplength, insertwidth, circleplotlength];

%% cumulative plot
axesinfo(4).position = [.13, .12 + 2*verticalseplength, 0.75*insertwidth, cumulplotlength];


%% cmapdata = [
%%     t1 t1 1
%%     1 1 t1
%%     1 t1 t1
%%     t2*[1 1 1]
%%            ];
%% %    t1 1 t1

t1 = 0.45;
t2 = 1;
cmapdata = [
     1 1 t1
     1 1 t1
    t1 t1 1
    t1 t1 1
           ];
%    t1 1 t1

colormap(cmapdata);

colors{1,1} = cmapdata(1,:);
colors{1,2} = cmapdata(2,:);
colors{2,1} = cmapdata(3,:);
colors{2,2} = cmapdata(4,:);

% automatically create postscript whenever
% figure is drawn

tmpdir = 'figmechturk_happy_wordchanges_comp001';
tmpfilename = sprintf('%s/%s_%s',picdir,tmpdir,file_tag);
tmpcommand = sprintf('mkdir -p %s',tmpdir);
%system(tmpcommand);

tmpfilenoname = sprintf('%s_noname',tmpfilename);
tmpfilenoname = strrep(tmpfilenoname,' ','_');
tmpfilenonamepdf = sprintf('%s.pdf',tmpfilenoname);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmpa1 = axes('position',axesinfo(1).position);

set(gcf,'DefaultAxesColor','none');
set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');

% tmpsym = {'ok-','sk-','dk-','vk-','^k-','>k-','<k-','pk-','hk-'};
% tmpsym = {'k-','r-','b-','m-','c-','g-','y-'};
% tmpsym = {'k-','k-.','k:','k--','r-','r-.','r:','r--'};
% tmplw = [ 1.5*ones(1,4), .5*ones(1,4)];

% main data goes here

load twitterwords_valence_vectors;
%load twitterwords;

%% cut vectors down
words = twitterwords(bigindex);
twitterwords_val_mean = twitterwords_val_mean(bigindex);
freq1 = freq1(bigindex);
freq2 = freq2(bigindex);

tmpwf1 = col(freq1/sum(freq1));
tmpwf2 = col(freq2/sum(freq2));

havg1 = sum(twitterwords_val_mean.*tmpwf1);
havg2 = sum(twitterwords_val_mean.*tmpwf2);

basemean = sum(twitterwords_val_mean.*tmpwf1);
deltas = (tmpwf2-tmpwf1).*(twitterwords_val_mean - basemean);
deltas = deltas/abs(sum(deltas))*100;

tmpmaxdelta = max(abs(min(deltas)),abs(max(deltas)));

%% summary circles
totalshifts = [0 0 ; 0 0];
%% entries correspond to position 
%% on figure according to normal
%% matrix indices; 
%% e.g., (1,1) is the upper left corner
%% and shows relatively positive words
%% being used less
for i=1:length(twitterwords_val_mean)
if ((twitterwords_val_mean(i)<basemean) & (deltas(i)>0))
    totalshifts(2,2) = totalshifts(2,2) + deltas(i);
elseif ((twitterwords_val_mean(i)>basemean) & (deltas(i)>0))
    totalshifts(1,2) = totalshifts(1,2) + deltas(i);
elseif ((twitterwords_val_mean(i)<basemean) & (deltas(i)<0))
    totalshifts(2,1) = totalshifts(2,1) + deltas(i);
else %% ((twitterwords_val_mean(i)>basemean) & (deltas(i)<0))
    totalshifts(1,1) = totalshifts(1,1) + deltas(i);
end
end

xlimmin = 0;
xlimmax = 0;

[tmp,ind] = sort(abs(deltas),'descend');

deltasort = deltas(ind(1:num));

tmp = deltas(ind(1:num));
deltamod = zeros(num,2);
tmpind = find(tmp>0);
deltamod(tmpind,1) = tmp(tmpind);
tmpind = find(tmp<0);
deltamod(tmpind,2) = tmp(tmpind);

[tmp2,ind2] = sort(abs(deltas),'descend');
tmp2 = deltas(ind2);
deltamod2 = zeros(length(deltas),2);
tmpind2 = find(tmp2>0);
deltamod2(tmpind2,1) = tmp2(tmpind2);
tmpind2 = find(tmp2<0);
deltamod2(tmpind2,2) = tmp2(tmpind2);

tmpbh = barh(deltamod,'stack');
set(gca,'ydir','reverse')

set(gca,'tickdir','out');

tmpch = get(tmpbh,'children');

set(tmpch{1},'facevertexcdata',1+zeros(num,1));
set(tmpch{2},'facevertexcdata',2+zeros(num,1));
set(tmpch{1},'cdatamapping','direct');
set(tmpch{2},'cdatamapping','direct');


ylim([0, (num+1)]);

hold on;

tmpycoords = 1:num;

%% output first 100 to terminal

for i=1:num
    % converts normalized x, y positions to appropriate values for axes
    %    [tmpXcoord,tmpYcoord] = normfigcoords(tmpxcoords(i),.005);
    tmpYcoord = tmpycoords(i);
    if (deltas(ind(i))<0)
        %        tmpYcoord = deltas(ind(i)) - 1;
        tmpXcoord = deltas(ind(i)) - 0.05*tmpmaxdelta;
    else
        %        tmpYcoord = deltas(ind(i)) + 1;
        tmpXcoord = deltas(ind(i)) + 0.05*tmpmaxdelta;
    end

    %% left to right, top to bottom
    %% 11, 12, 21, 22
    if ((twitterwords_val_mean(ind(i))>basemean) & (deltas(ind(i))<0))
        tmpstr = sprintf('+\\downarrow{%s}',words{ind(i)});
        %        tmpcolor = 'r';
        tmpcolor = [255 25 25]/255;
        tmpcolor = .05*[1 1 1];

        tmpfd = get(tmpch{2},'facevertexcdata');
        tmpfd(i)=1;
        set(tmpch{2},'facevertexcdata',tmpfd);

    elseif ((twitterwords_val_mean(ind(i))>basemean) & (deltas(ind(i))>0))
        tmpstr = sprintf('%s +\\uparrow',words{ind(i)});
        %        tmpcolor = 0.5*[1 1 1];
        tmpcolor = [255 25 25]/255;
        tmpcolor = .05*[1 1 1];

        tmpfd = get(tmpch{1},'facevertexcdata');
        tmpfd(i)=2;
        set(tmpch{1},'facevertexcdata',tmpfd);

    elseif ((twitterwords_val_mean(ind(i))<basemean) & (deltas(ind(i))<0))
        tmpstr = sprintf('-\\uparrow{%s}',words{ind(i)});
        %        tmpcolor = 'k';
        tmpcolor = [25 25 190]/255;
        tmpcolor = .5*[1 1 1];
    
        tmpfd = get(tmpch{2},'facevertexcdata');
        tmpfd(i)=3;
        set(tmpch{2},'facevertexcdata',tmpfd);
    elseif ((twitterwords_val_mean(ind(i))<basemean) & (deltas(ind(i))>0))
        tmpstr = sprintf('%s -\\downarrow',words{ind(i)});
        %        tmpcolor = 'b';
        tmpcolor = [25 25 190]/255;
        tmpcolor = .5*[1 1 1];
        
        tmpfd = get(tmpch{1},'facevertexcdata');
        tmpfd(i)=4;
        set(tmpch{1},'facevertexcdata',tmpfd);
    end

    fprintf(1,'%d: %s\n',i,tmpstr);
    
    tmpth(i) = text(tmpXcoord,tmpYcoord,tmpstr,...
                    'Fontsize',14, ...
                    'fontname','helvetica','color',tmpcolor);
    
        if (twitterwords_val_mean(ind(i))>basemean)
    %        set(tmpth(i),'color',[25 25 190]/255)
        else 
    %        set(tmpth(i),'color',[255 25 25]/255);
    %        set(tmpth(i),'fontangle','italic');
        end

    if (deltas(ind(i))<0)
        set(tmpth(i),'horizontalalignment','right');
    end
end

tmpxlim = get(gca,'xlim');
xlimmin = tmpxlim(1);
xlimmax = tmpxlim(2);

screenfactor = mainaxeswidth/(xlimmax-xlimmin);

%% relative position of zero

leftratios = zeros(size(deltasort));
rightratios = zeros(size(deltasort));

for i=1:num
    %% get extent of text
    extent = get(tmpth(i),'extent');
    
    %% screen width of text
    axestextwidths(i) = extent(3)*screenfactor;

    if (deltasort(i)<0)
        xlimminlocal = extent(1)+extent(3);
    else
        xlimmaxlocal = extent(1);
    end
end

indleft = find(deltasort<0);
indright = find(deltasort>=0);

alpha = 0.12;
tol = 0.0001;

lambda = 10;
lambdastep = 10;
[a,inda] = max([-deltasort(indleft)*screenfactor*lambda + axestextwidths(indleft)']);
[b,indb] = max([ deltasort(indright)*screenfactor*lambda + ...
                 axestextwidths(indright)']);
nsteps = 1;
while ((abs(a+b - mainaxeswidth*(1-alpha))>tol*mainaxeswidth) & (nsteps < 100))
    nsteps = nsteps + 1;
    lambdastep = lambdastep/2;
    if (a+b > mainaxeswidth*(1-alpha))
        lambda = lambda - lambdastep;
    else
        lambda = lambda + lambdastep;
    end
    [a,inda] = max([-deltasort(indleft)*screenfactor*lambda + axestextwidths(indleft)']);
    [b,indb] = max([ deltasort(indright)*screenfactor*lambda + ...
                     axestextwidths(indright)']);
end

set(gca,'fontsize',14);

% adjust scaling of limits
xlim(1/lambda*[xlimmin, xlimmax]);

%% now correct for left-right placement
tmpxlim = get(gca,'xlim');
xlimmin = tmpxlim(1);
xlimmax = tmpxlim(2);

screenfactor = mainaxeswidth/(xlimmax-xlimmin);

[a2,inda2] = max([-deltasort(indleft)*screenfactor + axestextwidths(indleft)']);
xlim([tmpxlim + (-a2/screenfactor - xlimmin - .5*alpha*mainaxeswidth/screenfactor)]);

tmpxlim = get(gca,'xlim');
set(gca,'xlim',[-1 1]*max(abs(tmpxlim)));

% fix up tickmarks
set(gca,'ytick',[1 5:5:num]);
% ridiculous (makes sure printing
% doesn't fiddle with xtick... thought i had solved this...
tmpxtick = get(gca,'xtick');
set(gca,'xtick',tmpxtick);

%% add circles

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
% tmplh = legend('','','');
% set(tmplh,'position',get(tmplh,'position')-[x y 0 0])
% change font
% tmplh = findobj(tmplh,'type','text');
% set(tmplh,'FontSize',18);
% remove box:
% legend boxoff

%% tmpxlab=xlabel('Per word happiness shift \delta{\ith}_{avg, \itr} (%)',...
tmpxlab=xlabel('Per word average happiness shift $\delta{h}_{{\rm avg}, r}$ (\%)',...
               'interpreter','latex',...
               'fontsize',14,'verticalalignment','top',...
               'horizontalalignment','center');

set(tmpxlab,'position',get(tmpxlab,'position')-[1,2, 0]);
%set(tmpxlab,'position',get(tmpxlab,'position')-[1, -0.5, 0]);
%% tmpylab=ylabel('Word rank \itr','fontsize',18,'verticalalignment','bottom');

tmpylab=ylabel('Word rank $r$','fontsize',18,'verticalalignment', ...
               'bottom','interpreter','latex','fontname','helvetica');

clear tmpstr;

tmpstr{1} = sprintf('$T_{\\rm ref}$  : %s ($h_{\\rm avg}$=%.2f)',...
                    title_text{1},havg1);
tmpstr{2} = sprintf('$T_{\\rm comp}$  : %s ($h_{\\rm avg}$=%.2f)',...
                    title_text{2},havg2);


[tmpXcoord,tmpYcoord] = normfigcoords(0.5,-.015);
tmpth = text('interpreter','latex',...
             'string',tmpstr,...
             'position',[tmpXcoord, tmpYcoord],...
             'Fontsize',15,...
             'verticalalignment','bottom',...
             'horizontalalignment','center',...
             'fontname','helvetica');



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% text size
%%%%%%%%%%%%%%%%%%%%%%%%%

tmpa2 = axes('position',axesinfo(2).position);

%% make into area:
totalmechturkwords = [sum(freq1); sum(freq2)].^.5;
maxN = 2.65*max(totalmechturkwords);

% left
d = totalmechturkwords(1)/maxN;
r = .5*d;
rvals(1) = r;
tmph = rectangle('position',[0.25-r, 0.25-r, d, d],'Curvature',[0,0]);
%% tmph = rectangle('position',[0.125, 0.25-r, 0.25, d],'Curvature',[0,0])
tmpcolor = 0.71*[1 1 1];
set(tmph,'facecolor',tmpcolor);
set(tmph,'edgecolor',0.31*[1 1 1]);

%% tmpstr = sprintf('+\\downarrow');
%% tmpXcoord = 0.0;
%% tmpYcoord = 0.95;
%% tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',16,'fontname','helvetica');

hold on;

% right
d = totalmechturkwords(2)/maxN;
r = .5*d;
rvals(2) = r;
tmph = rectangle('position',[0.75-r, 0.25-r, d, d],'Curvature',[0,0]);
%% tmph = rectangle('position',[0.625, 0.25-r, 0.25, d],'Curvature',[0,0]);
tmpcolor = 0.71*[1 1 1];
set(tmph,'facecolor',tmpcolor);
set(tmph,'edgecolor',0.31*[1 1 1]);

%% tmpstr = sprintf('+\\uparrow');
%% tmpXcoord = 1;
%% tmpYcoord = 0.95;
%% tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',16,'fontname', ...
%%             'helvetica','horizontalalignment','right');

%% dividing line

tmpx = 0:.01:1;
tmpxmid = mean([0.25+rvals(1); 0.75-rvals(2)]);
% tmph = plot(tmpx,0.5*ones(size(tmpx)),'k--');
% set(tmph,'color',0.5*[1 1 1]);
% hold on;

tmph = plot(tmpxmid*ones(size(tmpx)),tmpx,'k-');
set(tmph,'color',0.1*[1 1 1]);
hold on;

%% 

xlim([0 1]);
ylim([0 .5]);

%% text heads

if (havg1 > havg2)
    [tmpXcoord1,tmpYcoord1] = normfigcoords(0.25,0.95);
    [tmpXcoord2,tmpYcoord2] = normfigcoords(0.75,0.95);
else
    [tmpXcoord1,tmpYcoord1] = normfigcoords(0.75,0.95);
    [tmpXcoord2,tmpYcoord2] = normfigcoords(0.25,0.95);
end

[tmpXcoord1,tmpYcoord1] = normfigcoords(0.5,1.05);

tmpstr{1} = 'Text size:';
tmpstr{2} = '$T_{\rm ref}$ \ \  $T_{\rm comp}$';
tmpth1 = text('string',tmpstr,...
              'position',[tmpXcoord1, tmpYcoord1],...
              'interpreter','latex',...
              'fontsize',12,...
              'verticalalignment','bottom',...
              'horizontalalignment','center',...
              'fontname','helvetica');


%% [tmpXcoord1,tmpYcoord1] = normfigcoords(0.75,0.95);
%% [tmpXcoord2,tmpYcoord2] = normfigcoords(0.25,0.95);
%% 
%% tmpstr = '$T_{\rm ref}$';
%% tmpth2 = text('interpreter','latex',...
%%              'string',tmpstr,...
%%              'position',[tmpXcoord2, tmpYcoord2],...
%%              'Fontsize',16,...
%%              'verticalalignment','bottom',...
%%              'horizontalalignment','center',...
%%              'fontname','helvetica');
%% 
%% 
%% tmpstr = '$T_{\rm comp}$';
%% tmpth1 = text('interpreter','latex',...
%%              'string',tmpstr,...
%%              'position',[tmpXcoord1, tmpYcoord1],...
%%              'Fontsize',16,...
%%              'verticalalignment','bottom',...
%%              'horizontalalignment','center',...
%%              'fontname','helvetica');

set(gca,'visible','off');


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% circle contribution plot
%%%%%%%%%%%%%%%%%%%%%%%%%

tmpa2 = axes('position',axesinfo(3).position);

%% make into area:
totalshifts_orig = totalshifts;

totalshifts = (abs(totalshifts)).^.5;

maxts = 2.65*max(abs(totalshifts(:)));

tmpx = 0:.01:1;
% tmph = plot(tmpx,0.5*ones(size(tmpx)),'k--');
% set(tmph,'color',0.5*[1 1 1]);
% hold on;

tmph = plot(0.5*ones(size(tmpx)),tmpx,'k-');
set(tmph,'color',0.1*[1 1 1]);
hold on;

% top left
d = abs(totalshifts(1,1))/maxts;
r = .5*d;
tmph = rectangle('position',[0.25-r, 0.75-r, d, d],'Curvature',[1,1]);
tmpcolor = colors{1,1};
set(tmph,'facecolor',tmpcolor);
set(tmph,'edgecolor',0.31*[1 1 1]);

tmpstr = sprintf('+\\downarrow');
tmpXcoord = 0.0;
tmpYcoord = 0.95;
tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',12,'fontname','helvetica');

hold on;

% top right
d = abs(totalshifts(1,2))/maxts;
r = .5*d;
tmph = rectangle('position',[0.75-r, 0.75-r, d, d],'Curvature',[1,1]);
tmpcolor = colors{1,2};
set(tmph,'facecolor',tmpcolor);
set(tmph,'edgecolor',0.31*[1 1 1]);

tmpstr = sprintf('+\\uparrow');
tmpXcoord = 1;
tmpYcoord = 0.95;
tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',12,'fontname', ...
             'helvetica','horizontalalignment','right');

% bottom left
d = abs(totalshifts(2,1))/maxts;
r = .5*d;
tmph = rectangle('position',[0.25-r, 0.25-r, d, d],'Curvature',[1,1]);
tmpcolor = colors{2,1};
set(tmph,'facecolor',tmpcolor);
set(tmph,'edgecolor',0.31*[1 1 1]);

tmpstr = sprintf('-\\uparrow');
tmpXcoord = 0.00;
tmpYcoord = 0.05;
tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',12,'fontname', ...
             'helvetica','horizontalalignment','left');
set(tmpth,'color',.5*[1 1 1]);

% bottom right
d = abs(totalshifts(2,2))/maxts;
r = .5*d;
tmph = rectangle('position',[0.75-r, 0.25-r, d, d],'Curvature',[1,1]);
tmpcolor = colors{2,2};
set(tmph,'facecolor',tmpcolor);
set(tmph,'edgecolor',0.31*[1 1 1]);

tmpstr = sprintf('-\\downarrow');
tmpXcoord = 1;
tmpYcoord = 0.05;
tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',12,'fontname', ...
             'helvetica','horizontalalignment','right');
set(tmpth,'color',.5*[1 1 1]);

xlim([0 1]);
ylim([0 1]);

set(gca,'visible','off');

[tmpXcoord1,tmpYcoord1] = normfigcoords(0.5,1.3);
clear tmpstr;
tmpstr = 'Balance:';
tmpth1 = text('string',tmpstr,...
              'position',[tmpXcoord1, tmpYcoord1],...
              'interpreter','latex',...
              'fontsize',12,...
              'verticalalignment','bottom',...
              'horizontalalignment','center',...
              'fontname','helvetica');

[tmpXcoord1,tmpYcoord1] = normfigcoords(0.5,1.05);
tmpstr = sprintf('%.0f : +%.0f',sum(totalshifts_orig(:,1)),sum(totalshifts_orig(:,2)));
tmpth1 = text('string',tmpstr,...
              'position',[tmpXcoord1, tmpYcoord1],...
              'fontsize',12,...
              'verticalalignment','bottom',...
              'horizontalalignment','center',...
              'fontname','helvetica');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% cumulative circle contribution plot
%%%%%%%%%%%%%%%%%%%%%%%%%

tmpa2 = axes('position',axesinfo(4).position);

tmpcumsum = cumsum(sum(deltamod2(1:end,:),2));
tmph = semilogy(tmpcumsum,(1:size(deltamod2,1)),'k.-');
set(tmph,'linewidth',0.5);
set(tmph,'markersize',5);

%% ylim([0, length(words)+1]);
ylim([0, 10^4]);
xlim(1.03*[min([tmpcumsum; 0]), max([tmpcumsum; 0])]);

set(gca,'linewidth',0.5);

min(tmpcumsum)

tmpxlim = get(tmpa2,'xlim');

hold on;
tmpx = linspace(tmpxlim(1),tmpxlim(2),1000);
plot(tmpx,num*ones(size(tmpx)),'k-');

%addlabel(labelin,-1,2.425,18);


set(gca,'fontsize',10);
set(gca,'ydir','reverse')

set(gca,'ytick',[1 10 100 1000 10000]);
set(gca,'xtick',[-100 0 100]);

%% set(gca,'yaxislocation','right');


tmpxlab=xlabel('$\sum_{i=1}^{r} \delta{h}_{{\rm avg}, i}$',...
               'interpreter','latex',...
               'fontsize',12,'verticalalignment','top');
%% tmpxlab=xlabel('\sffamily$\Sigma_{{i}=1}^{r} \delta{h}_{{\rm avg}, i}$ (\%)',...
%%               'interpreter','latex',...
%%               'fontsize',16,'verticalalignment','top');
%% set(tmpxlab,'position',get(tmpxlab,'position')-[0, -4, 0]);

% set(tmpxlab,'position',get(tmpxlab,'position') - [0 .1 0]);
% set(tmpylab,'position',get(tmpylab,'position') + [2 10 0]);

%if (sum(deltas)<0)
%  [tmpXcoord,tmpYcoord] = normfigcoords(.5,.8);
%else
%  [tmpXcoord,tmpYcoord] = normfigcoords(.5,.2);
%end
%

[tmpXcoord,tmpYcoord] = normfigcoords(1.8,.02);

clear tmpstr
tmpstr{1} = sprintf('Mean happiness:');
tmpstr{2} = sprintf('%4.2f vs. %4.2f',sum(tmpwf2.*twitterwords_val_mean), ...
                    sum(tmpwf1.*twitterwords_val_mean));
%% disp(tmpstr);
%%
%% tmpth = text(tmpXcoord,tmpYcoord,tmpstr,'Fontsize',16,'fontname','helvetica','rotation',90);


%%%%%%%%%%%%%%
%%% print out

% automatic creation of postscript
% without name/date
%psprintcpdf(tmpfilenoname);

%% tmpcommand = sprintf('open -a preview');
%% system(tmpcommand);

 %% tmpcommand = sprintf('open %s',tmpfilenonamepdf);
 %% system(tmpcommand);

 % title('','fontsize',24)
 % 'horizontalalignment','left');
 % tmpxl = xlabel('','fontsize',24,'verticalalignment','top');
 % set(tmpxl,'position',get(tmpxl,'position') - [ 0 .1 0]);
 % tmpyl = ylabel('','fontsize',24,'verticalalignment','bottom');
 % set(tmpyl,'position',get(tmpyl,'position') - [ 0.1 0 0]);
 % title('','fontsize',24)

 % converts normalized x, y positions to appropriate values for axes
 % [tmpXcoord,tmpYcoord] = normfigcoords(,);
 % text(tmpXcoord,tmpYcoord,'','Fontsize',14);

 % name label
 %% tmpt = pwd;
 %% tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
 %% 
 %% [tmpXcoord,tmpYcoord] = normfigcoords(1.05,.05);
 %% text(tmpXcoord,tmpYcoord,tmpnamememo,'Fontsize',6,'rotation',90);
 %% 
 %% [tmpXcoord,tmpYcoord] = normfigcoords(1.1,.05);
 %% datenamer(tmpXcoord,tmpYcoord,90);
 %% 
 % [tmpXcoord,tmpYcoord] = normfigcoords(.5,.05);
 % datename(tmpXcoord,tmpYcoord);

 % [tmpXcoord,tmpYcoord] = normfigcoords(.5,.05);
 % datename2(tmpXcoord,tmpYcoord); % 2 rows

 % automatic creation of postscript
 % psprint(tmpfilename);

clear tmp*

totalshifts;
totalshifts_orig;

%more on;
