% Copyright : This code is written by Gerasimos Arvanitis {arvanitis@ece.upatras.gr}
%              Electrical and Computer Engineering, University of Patras. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion of this copyright information.

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
