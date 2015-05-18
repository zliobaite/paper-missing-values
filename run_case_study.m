%This code is to support the article:

%Zliobaite, I., Hollmen, J. and Junninen, H. (2014).
%Predictive models tolerant to massively missing data: a case study in solar radiation prediction. Currently under review at Atmospheric Environment, Elsevier.

%The data and the code can be used for research purposes, provided that the above article is cited.

%This code is available from http://users.ics.aalto.fi/indre/smear.zip

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
%ertype2 = 'ianew';
components = 18;
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
labels = labels*100; %in %

%keep only times with daylight
indkeep = labels>0; 
data = data(indkeep,:);
thrad = thrad(indkeep);
labels = labels(indkeep);
dates = dates(indkeep,:);

[n,k] = size(data);
disp(['data size ',num2str(n),' x ',num2str(k)]);

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

data2(isnan(data2))=0;

ss = cov(data1);
    
%[comp,sc,lat] = pca(data1);
%[[1:36]' cumsum(lat/sum(lat))]
comp = pca_reg(data1);

variable_contributions = zeros(k,1);
%coreelation based
for sk2 = 1:k
	R=corrcoef(data1(:,sk2),y1);
    variable_contributions(sk2) = R(2,1);
end;
[i,var_sorted] = sort(abs(variable_contributions),'descend');

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



%all_models = [betaALL betaALLreg betaFSE betaFSEreg betaPCA betaPCAreg betaPLS];
    
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

errors_train_mse = [];
errors_train_R2 = [];
errors_test = [];
errors_testA = [];
errors_testM = [];
%ia_train = [];
%ia_test = [];
%ia_testA = [];
%ia_testM = [];


y2_orig = standardize_back(y2,ymean,ystd);
y1_orig = standardize_back(y1,ymean,ystd);


for sk3=1:8
    
    errors_train_mse = [errors_train_mse error_reg(y1_orig,predictions1(:,sk3),ertype)];
    errors_train_R2 = [errors_train_R2 error_reg(y1_orig,predictions1(:,sk3),'R2')];
    
    errors_test = [errors_test error_reg(y2_orig,predictions2(:,sk3),ertype)];
    errors_testA = [errors_testA error_reg(y2_orig(ind2),predictions2(ind2,sk3),ertype)];
    errors_testM = [errors_testM error_reg(y2_orig(indnan2),predictions2(indnan2,sk3),ertype) ];
    
    %ia_train = [ia_train error_reg(y1_orig,predictions1(:,sk3),ertype2)];
    %ia_test = [ia_test error_reg(y2_orig,predictions2(:,sk3),ertype2)];
    %ia_testA = [ia_testA error_reg(y2_orig(ind2),predictions2(ind2,sk3),ertype2)];
    %ia_testM = [ia_testM error_reg(y2_orig(indnan2),predictions2(indnan2,sk3),ertype2)];
end;

residuals = abs(predictions2 - y2_orig*ones(1,9));

betabeta = [betaALL'*betaALL betaALLreg'*betaALLreg betaFSE'*betaFSE betaFSEreg'*betaFSEreg betaPCA'*betaPCA betaPCAreg'*betaPCAreg betaPLS'*betaPLS];
 
DDnow = [];
DDnow = [DDnow betaALL'*(ss-eye(k))*betaALL];
DDnow = [DDnow betaALLreg'*(ss-eye(k))*betaALLreg];
DDnow = [DDnow betaFSE'*(ss-eye(k))*betaFSE];
DDnow = [DDnow betaFSEreg'*(ss-eye(k))*betaFSEreg];
DDnow = [DDnow betaPCA'*(ss-eye(k))*betaPCA];
DDnow = [DDnow betaPCAreg'*(ss-eye(k))*betaPCAreg];
DDnow = [DDnow betaPLS'*(ss-eye(k))*betaPLS];
DDnow = -DDnow;

DDnow = DDnow*ystd*ystd;

varyy = [betaALL'*ss*betaALL betaALLreg'*ss*betaALLreg betaFSE'*ss*betaFSE betaFSEreg'*ss*betaFSEreg betaPCA'*ss*betaPCA betaPCAreg'*ss*betaPCAreg betaPLS'*ss*betaPLS];


record_ci = [[1:36]' betaALL betaALLreg betaFSE betaFSEreg betaPCA betaPCAreg betaPLS ones(k,1)*sqrt(errors_train_mse*n/(n-k))*1.645];
record_ci(:,2:8) = record_ci(:,2:8)*ystd;
for sk=1:7
    ind = find(record_ci(:,sk+1)==0);
    if ~isempty(ind)
        record_ci(ind,sk+8)=0;
    end;
end;
csvwrite('conf_int.csv',record_ci);


%cross-validation
nn = size(data1,1);
folds = 10;
rng(1); %random seed
indices = crossvalind('Kfold', nn, folds);

predALL = zeros(nn,1);
predALLreg = zeros(nn,1);
predFSE = zeros(nn,1);
predFSEreg = zeros(nn,1);
predPCA = zeros(nn,1);
predPCAreg = zeros(nn,1);
predPLS = zeros(nn,1);

for sk=1:folds
    
    test_range = find(indices==sk);
    train_range = intersect([1:nn],test_range);
    
    data_train = data1(train_range,:);
    data_test = data1(test_range,:);
    labels_train = y1(train_range);
    
    %[comp,sc,lat] = pca(data_train);
    comp = pca_reg(data_train);

    variable_contributions = zeros(k,1);
    %coreelation based
    for sk2 = 1:k
        R=corrcoef(data_train(:,sk2),labels_train);
        variable_contributions(sk2) = R(2,1);
    end;
    [i,var_sorted] = sort(abs(variable_contributions),'descend');
    

    betaALL = reg_regression_train(data_train,labels_train,0);
    
    betaALLreg = reg_regression_train(data_train,labels_train,200);

    modelFSE = reg_regression_train(data_train(:,var_sorted(1:components)),labels_train,0);
    betaFSE = zeros(k,1);
    betaFSE(var_sorted(1:components)) = modelFSE;
    
    modelFSEreg = reg_regression_train(data_train(:,var_sorted(1:components)),labels_train,200);
    betaFSEreg = zeros(k,1);
    betaFSEreg(var_sorted(1:components)) = modelFSEreg;
    
    modelPCA = reg_regression_train(data_train*comp(:,1:components),labels_train,0);
    betaPCA = comp(:,1:components)*modelPCA;

    modelPCAreg = reg_regression_train(data_train*comp(:,1:components),labels_train,200);
    betaPCAreg = comp(:,1:components)*modelPCAreg;

    [W,P,Q] = nipals_train_batch_nomean(data_train,labels_train,components);
    betaPLS = W*inv(P'*W)*Q;
    
    predALL(test_range) = data_test*betaALL;
    predALLreg(test_range) = data_test*betaALLreg;
    predFSE(test_range) = data_test*betaFSE;
    predFSEreg(test_range) = data_test*betaFSEreg;
    predPCA(test_range) = data_test*betaPCA;
    predPCAreg(test_range) = data_test*betaPCAreg;
    predPLS(test_range) = data_test*betaPLS;
    
    predALL(test_range) = standardize_back(predALL(test_range),ymean,ystd);
    predALLreg(test_range) = standardize_back(predALLreg(test_range),ymean,ystd);
    predFSE(test_range) = standardize_back(predFSE(test_range),ymean,ystd);
    predFSEreg(test_range) = standardize_back(predFSEreg(test_range),ymean,ystd);
    predPCA(test_range) = standardize_back(predPCA(test_range),ymean,ystd);
    predPCAreg(test_range) = standardize_back(predPCAreg(test_range),ymean,ystd);
    predPLS(test_range) = standardize_back(predPLS(test_range),ymean,ystd);
end;

errors_cv = [error_reg(y1_orig,predALL,ertype) error_reg(y1_orig,predALLreg,ertype) error_reg(y1_orig,predFSE,ertype) error_reg(y1_orig,predFSEreg,ertype) error_reg(y1_orig,predPCA,ertype) error_reg(y1_orig,predPCAreg,ertype) error_reg(y1_orig,predPLS,ertype) error_reg(y1_orig,predictions1(:,8),ertype)];
%ia_cv = [error_reg(y1_orig,predALL,ertype2) error_reg(y1_orig,predALLreg,ertype2) error_reg(y1_orig,predFSE,ertype2) error_reg(y1_orig,predFSEreg,ertype2) error_reg(y1_orig,predPCA,ertype2) error_reg(y1_orig,predPCAreg,ertype2) error_reg(y1_orig,predPLS,ertype2) error_reg(y1_orig,predictions1(:,8),ertype2)];
    

MM = [errors_test; errors_testA; errors_testM; errors_cv; errors_train_mse; errors_train_R2];
MMIA = [ia_test; ia_testA; ia_testM; ia_cv];

disp('Table 3');
disp('10-fold cross validation errors (RMSE) on the training dataset and deterioration index (d)');
disp( '         ALL         rALL          SEL         rSEL          PCA         rPCA          PLS');
disp(['RMSE    ',num2str(round(MM(4,1:end-1)*10)/10)]);
disp(['d index ',num2str(round(DDnow))]);
disp('--');

disp('Table 5');
disp('Prediction errors (RMSE) on the testing dataset (%)');
disp( '             ALL         rALL          SEL         rSEL          PCA         rPCA          PLS          NAI');
disp(['full set    ',num2str(round(MM(1,:)*10)/10)]);
disp(['non-missing ',num2str(round(MM(2,:)*10)/10)]);
disp(['missing     ',num2str(round(MM(3,:)*10)/10)]);
disp('--');

% disp('Table 13');
% disp('10-fold cross validation results (IA) on the training dataset');
% disp( '         ALL         rALL          SEL         rSEL          PCA         rPCA          PLS');
% disp(['IA    ',num2str(round(MMIA(4,:)*100)/100)]);
% disp('--');
% disp('Table 14');
% disp('Prediction results (IA) on the testing dataset');
% disp( '             ALL         rALL          SEL         rSEL          PCA         rPCA          PLS          NAI');
% disp(['full set    ',num2str(round(MMIA(1,:)*100)/100)]);
% disp(['non-missing ',num2str(round(MMIA(2,:)*100)/100)]);
% disp(['missing     ',num2str(round(MMIA(3,:)*100)/100)]);

disp('Table B1');
disp('Fitness statistics of the regression models on the training data');
disp( '         ALL         rALL          SEL         rSEL          PCA         rPCA          PLS    NAI');
disp(['RMSE    ',num2str(round(MM(5,:)*10)/10)]);
disp(['R2 ',num2str(round(MM(6,:)*1000)/1000)]);

%histograms of residuals
%Figure 5
v = 1;
for sk=1:12
    v = [v; v(end)*2];
end;
v = [-1;v];

histos = zeros(length(v)-2,8);
for sk=1:9
    h = histc(residuals(:,sk),v);
    h(2) = h(2) + h(1);
    h = h(2:end-1);
    sum(h)
    histos(:,sk) = h/sum(h);
end;

disp('Figure 5');
disp('Analysis of residuals.');
disp( 'h         ALL         rALL          SEL         rSEL          PCA         rPCA          PLS    NAI');
disp([v(3:end)/2 round(histos*100)]);