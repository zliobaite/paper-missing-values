%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Last updated: 2013 09 21
%---------------------------------

function [W,P,Q] = nipals_train_batch_nomean(data,labels,number_of_components)

[n,p] = size(data);
X = data;
y = labels;


W = [];
P = [];
Q = [];

for sk=1:number_of_components
    
    w = X'*y;
    w = w/sqrt(w'*w); %optional 
    
    t = X*w;

    q = t'*y/(t'*t);
    p = X'*t/(t'*t);
        

    X = X - t*p'; %deflation
    y = y - t*q;
    
    W = [W w];
    P = [P p];
    Q = [Q; q];
    
end;

    





