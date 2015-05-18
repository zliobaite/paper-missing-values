%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Last updated: 2013 09 21
%---------------------------------

function [data,labels] = remove_missing_values(data,labels)

[n,p] = size(data);
keep = [1:n];
for sk=1:p
    ind = find(isnan(data(:,sk))==0);
    keep = intersect(keep,ind);
end;

data = data(keep,:);
labels = labels(keep);
labels = labels(:);