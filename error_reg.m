%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Last updated: 2014 01 28 (added Index of Agreement)
%---------------------------------

function err = error_reg(ytrue,ypredicted,ertype)

switch ertype
    case 'mse'
        err = mean( (ytrue - ypredicted).^2 ) ;
    case 'rmse'
        err = sqrt(mean( (ytrue - ypredicted).^2 ) );
    case 'mae'
        err = mean(abs(ytrue-ypredicted));
    case 'ia'
        top = sum(abs(ytrue - ypredicted).^2);
        mytrue = mean(ytrue);
        bot = sum((abs(ypredicted - mytrue) + abs(ytrue - mytrue)).^2);
        err = 1 - (top/bot);
    case 'ianew'
        top = sum(abs(ytrue - ypredicted));
        mytrue = mean(ytrue);
        bot = sum(2*abs(ytrue - mytrue));
        err = 1 - (top/bot);
    case 'R2'
        top = sum((ytrue - ypredicted).^2);
        mytrue = mean(ytrue);
        bot = sum((ytrue - mytrue).^2);
        err = 1 - (top/bot);
    otherwise
        disp('no error assigned');
end;