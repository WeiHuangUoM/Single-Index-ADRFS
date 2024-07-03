pg = 1:3;
Kg2 = 1:3;

[pg1,pg2,Kg2] = ndgrid(pg,pg,Kg2);
pg1 = pg1(:);
pg2 = pg2(:);
Kg2 = Kg2(:);
lK = length(pg1);

CV_all = zeros(1,lK);



for j = 1:lK
    p1 = pg1(j);
    p2 = pg2(j);
    K2 = Kg2(j);

    filename=sprintf('CV-tc0-%d-%d-%d-30-May-2024.mat', p1,p2,K2);
    try  load(filename)
        %load('/Users/zyl/Library/CloudStorage/OneDrive-Personal/Non_CTE20220226 _2/Matlab Code/real-data/filename')
        
        CV_all(j) = CV;
    catch ME
        filename=sprintf('CV-tc0-%d-%d-%d-01-Jun-2024.mat', p1,p2,K2);
        load(filename)
        CV_all(j) = CV;
    end
end

    index = find(CV_all == min(CV_all),1);
    p1 = pg1(index);
    p2 = pg2(index);
    K2 = Kg2(index);
