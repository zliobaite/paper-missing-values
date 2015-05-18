%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Last updated: 2013 09 21
%---------------------------------
ssdata = load('data_smearii.csv');
thrad = load('data_theoretical_radiation.csv');

dates = ssdata(:,1:6);
data = ssdata(:,[7:end]);

[n,p] = size(data);

missing = zeros(p,1);
durations = zeros(p,1);
std_durations = zeros(p,1);
for sk=1:p
    data_now = data(:,sk);
    inn = isnan(data_now);
    missing(sk) = sum(inn)/n;
    if ((sum(inn)>0) & (sum(inn)~=n))
        ind1 = find(diff(inn)==1);
        ind2 = find(diff(inn)==-1);
        if ind1(1)>ind2(1)
            ind1 = [0 ; ind1];
        end;
        if length(ind1)>length(ind2)
            ind2 = [ind2 ; n];
        end;
        durations(sk) = mean(ind2 - ind1);
        std_durations(sk) = std(ind2 - ind1);
    end;
end;

dd = isnan(data); sdd = sum(dd,2);
dd = dd(:);
disp(['Percentage of missing values ',num2str(mean(dd)*100),' %']);
disp('--');
disp(['Figure 2']);
disp('(a)');

ddt = sum(isnan(data(:,1:end-1)),2);
N=histc(ddt,[0:36]); %hist
N = [N(1:10); sum(N(11:end))];
disp(['number of missing sensors       number of observations']);
for sk=0:9
    disp([num2str(sk),'  ',num2str(round(N(sk+1)*1000/sum(N))/10),'%']);
end;
disp(['10+  ',num2str(round(N(11)*1000/sum(N))/10),'%']);



missing_share = zeros(p-1,5);
missing_share(:,1) = [0:p-2];
   
%prepare labels
labels = ssdata(:,end); labels(labels<0)=0;
ssdata = ssdata(:,7:end-1);

labels(labels<0)=0; %removes negative observations
labels = labels ./ thrad(:,end);

labels(isnan(labels))=0;
labels(labels==Inf)=0;
labels(labels>1)=1; 
labels(labels<0)=0;

%only light days
indkeep = labels>0; 
%indkeep = [1:length(thrad)];
data = ssdata(indkeep,:);
labels = labels(indkeep);
dates = dates(indkeep,:);

[N,P] = size(ssdata);
[n,p] = size(data);

no = sum(isnan(ssdata),1);
[i,j] = sort(no);



for sk=1:p
    mm = sum(isnan(ssdata(:,j(1:p-sk+1))),2);
    missing_share(sk,2) = length(find(mm>0))/N;
    
    ind1 = find(isnan(labels)==0);
    ind2 = find(isnan(data(:,j(p-sk+1)))==0);
    ind = intersect(ind1,ind2);
    
    missing_share(sk,5) = abs(corr(labels(ind),data(ind,j(p-sk+1))));
    missing_share(sk,4) = 1- length(ind2)/n;
    
end;   
missing_share(:,3) = 1 - missing_share(:,2);


disp('(b)');
disp(['number of remaining sensors       % of complete data']);
for sk=1:p
    disp([num2str(p-sk+1),'  ',num2str(round(missing_share(sk,3)*1000)/10),'%']);
end;

disp('(c)');
disp('dark times removed');
disp(['missing rate %     absolute correlation with target']);
for sk=1:p
    disp([num2str(round(missing_share(sk,4)*1000)/10),'%  ',num2str(round(missing_share(sk,5)*1000)/1000)]);
end;

disp('Figure 4');
disp(['number of input sensors       mean duration (h)      standard veviation of the duration']);
for sk=1:p+1
    disp([num2str(sk),'  ',num2str(round(durations(sk)*10/2)/10),' ',num2str(round(std_durations(sk)*10/2)/10)]);
end;
%dividion by 2 because data is in 30 min intervals
