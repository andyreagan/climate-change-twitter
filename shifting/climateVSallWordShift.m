load allDates29.mat
load climateWFM2.mat
load completeDailyWFM.mat
load twitterwords_valence_vectors.mat

close all;

%parameters (if only one day, put in same day twice)
firstDay = '2008-12-28';
lastDay = '2008-12-28';

%data
dateList = allDates29;
WFM1 = climateWFM2;
WFM2 = completeDailyWFM_cut';
saveeps = 1;

%get rid of neutral words
ishappy_indices = find((twitterwords_val_mean > 0) &...
    ((twitterwords_val_mean >= 6) | (twitterwords_val_mean <= 4)));

%get rid of nigga and niggas
nigga = find(ismember(twitterwords,'nigga'));
niggaIdx = find(ishappy_indices==nigga);
ishappy_indices(niggaIdx)=[];

niggas = find(ismember(twitterwords,'niggas'));
niggasIdx = find(ishappy_indices==niggas);
ishappy_indices(niggasIdx)=[];

indices = ishappy_indices;

%get indices of given date range
dateIDXs = find(ismember(allDates29,{firstDay, lastDay}));

%if only one day is given
if length(dateIDXs) == 1
    vec1 = nansum(WFM1(:,dateIDXs),2);
    vec2 = nansum(WFM2(:,dateIDXs),2);
    datetitle = firstDay;
else
    %if multiple days are given
    vec1 = nansum(WFM1(:,dateIDXs(1):dateIDXs(2)),2);
    vec2 = nansum(WFM2(:,dateIDXs(1):dateIDXs(2)),2);
    datetitle = strcat(firstDay,',',lastDay);
end

%word shift
subplot(2,1,1)
figmechturk_happy_wordchanges_comp001(indices,vec2,vec1,'name1', ...
    {firstDay,'Climate'},50)

if saveeps
    %save as firstDay.pdf
    psprintcpdf(firstDay)
end
    