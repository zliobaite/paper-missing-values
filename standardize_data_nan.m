%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Last updated: 2013 09 21
%---------------------------------

function new_data = standardize_data_nan(data,data_mean,data_std)
    
    [n,k] = size(data);
    
     new_data = (data - ones(n,1)*data_mean) ./ (ones(n,1)*data_std);