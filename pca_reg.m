%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Last updated: 2013 09 22

function R = pca_reg(data)

[R,~] = eig(cov(data));

R = fliplr(R);