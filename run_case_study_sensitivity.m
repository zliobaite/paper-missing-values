%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2013).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/data_code_smear.zip

%Mailto: zliobaite@gmail.com 
%Updated: 2014 01 28 (added Index of Agreement)
%Last updated: 2014 06 16 (reverted to RMSE, removed Index of Agreement)
%---------------------------------

ssdata = load('data_smearii.csv');
thrad = load('data_theoretical_radiation.csv');

%remove missing target variables
indkeep = find(isnan(ssdata(:,end))==0);
ssdata = ssdata(indkeep,:);
thrad = thrad(indkeep,end);

sunshine = zeros(size(thrad));
sunshine(thrad>0) = 1;

ertype = 'rmse';
components = 10;
docap = 0; %cap predictions

%prepare data
dates = ssdata(:,1:6);
data = ssdata(:,7:end-1); 

%prepare labels
labels = ssdata(:,end); labels(labels<0)=0;
labels(labels<0)=0; %removes negative observations
labels = labels ./ thrad;
labels(isnan(labels))=0;
labels(labels==Inf)=0;
labels(labels>1)=1; 
labels(labels<0)=0;
labels = labels * 100; %in %

%only light days
indkeep = labels>0; 
data = data(indkeep,:);
thrad = thrad(indkeep);
labels = labels(indkeep);
dates = dates(indkeep,:);

[n,k] = size(data);

ind = intersect( intersect( find(dates(:,1)==2008) , find(dates(:,2)==4) ) , find(dates(:,3)==15) );
train_range = [1:ind(1)-1];
test_range = [ind(1):n-2];

data1 = data(train_range,:);
y1 = labels(train_range);

data2 = data(test_range,:);
y2 = labels(test_range);   


[data1,y1] = remove_missing_values(data1,y1);


[dmean,dstd,data1]= standardize_data_nan_train(data1);
data2= standardize_data_nan(data2,dmean,dstd);
[ymean,ystd,y1]= standardize_data_nan_train(y1);
y2= standardize_data_nan(y2,ymean,ystd);

indnan2 = sum(isnan(data2),2);
indnan2 = find(indnan2>=1);
ind2 = setdiff([1:size(data2,1)],indnan2);

data2(isnan(data2))=0; %replacing missing values

ss = cov(data1);
%[comp,sc,lat] = pca(data1);
comp = pca_reg(data1);

variable_contributions = zeros(k,1);
%coreelation based
for sk2 = 1:k
	R=corrcoef(data1(:,sk2),y1);
    variable_contributions(sk2) = R(2,1);
end;
[i,var_sorted] = sort(abs(variable_contributions),'descend');
   

ET = [];
%IAT = [];

%sensitivity analysis
for components=1:k

    betaALL = reg_regression_train(data1,y1,0);

    betaALLreg = reg_regression_train(data1,y1,200);

    modelFSE = reg_regression_train(data1(:,var_sorted(1:components)),y1,0);
    betaFSE = zeros(k,1);
    betaFSE(var_sorted(1:components)) = modelFSE;

    modelFSEreg = reg_regression_train(data1(:,var_sorted(1:components)),y1,200);
    betaFSEreg = zeros(k,1);
    betaFSEreg(var_sorted(1:components)) = modelFSEreg;

    modelPCA = reg_regression_train(data1*comp(:,1:components),y1,0);
    betaPCA = comp(:,1:components)*modelPCA;

    modelPCAreg = reg_regression_train(data1*comp(:,1:components),y1,200);
    betaPCAreg = comp(:,1:components)*modelPCAreg;

    [W,P,Q] = nipals_train_batch_nomean(data1,y1,components);
    betaPLS = W*inv(P'*W)*Q;

    predictions1 = [data1*betaALL data1*betaALLreg data1*betaFSE data1*betaFSEreg data1*betaPCA data1*betaPCAreg data1*betaPLS ones(size(y1))*mean(y1) [0;0; y1(1:end-2)]]; 
    predictions2 = [data2*betaALL data2*betaALLreg data2*betaFSE data2*betaFSEreg data2*betaPCA data2*betaPCAreg data2*betaPLS ones(size(y2))*mean(y2) [y1(end-1:end); y2(1:end-2)]];

    %standardise back
    predictions1 = standardize_back(predictions1,ymean,ystd);
    predictions2 = standardize_back(predictions2,ymean,ystd);



    if docap
        for sk3=1:7
            ind = predictions1(:,sk3)>1;
            predictions1(ind,sk3)=1;
            ind = predictions1(:,sk3)<0;
            predictions1(ind,sk3)=0;

            ind = predictions2(:,sk3)>1;
            predictions2(ind,sk3)=1;
            ind = predictions2(:,sk3)<0;
            predictions2(ind,sk3)=0;

        end;
    end;

    errors_train = [];
    errors_test = [];
    errors_testA = [];
    errors_testM = [];
    %ia_test = [];

    y2_orig = standardize_back(y2,ymean,ystd);
    y1_orig = standardize_back(y1,ymean,ystd);
    for sk3=1:8

        errors_train = [errors_train error_reg(y1_orig,predictions1(:,sk3),ertype)];
        errors_test = [errors_test error_reg(y2_orig,predictions2(:,sk3),ertype)];
        %ia_test = [ia_test error_reg(y2_orig,predictions2(:,sk3),'ianew')];
        
        errors_testA = [errors_testA error_reg(y2_orig(ind2),predictions2(ind2,sk3),ertype)];
        errors_testM = [errors_testM error_reg(y2_orig(indnan2),predictions2(indnan2,sk3),ertype) ];
    end;


ET = [ET; [components errors_test]];
%IAT = [IAT; [components ia_test]];
end;

disp('Figure 5');
disp('Prediction error as a function of components preserved');

disp('k     ALL         rALL          SEL         rSEL          PCA         rPCA          PLS          NAI');
for sk=1:36
    disp([num2str(sk),'   ',num2str(round(ET(sk,2:end)*10)/10)]);
end;
disp('k - components preserved');


% disp('Figure 5');
% disp('Prediction error  (in IA) as a function of components preserved');
% 
% disp('k     ALL         rALL          SEL         rSEL          PCA         rPCA          PLS          NAI');
% for sk=1:36
%     disp([num2str(sk),'   ',num2str(round(IAT(sk,2:end)*1000)/1000)]);
% end;
% disp('k - components preserved');