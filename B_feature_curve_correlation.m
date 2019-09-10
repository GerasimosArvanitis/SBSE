clear all

n = zeros(15,20);

for i = 1:15
    newfile = strcat(int2str(i),'_coefficient.csv');
    n(i,:) = xlsread(newfile);
end

for i = 1:15
    for j = 1:15
        corr(i,j) = corr2(n(i,:), n(j,:));
    end
end
       

norm_corr = 1-(corr - min(min(corr)))./ ( max(max(corr)) - min(min(corr)) );

csvwrite('similarities_of_feature_curves.csv',norm_corr);