function  [ECdist] = gngSessionSum2_PCA_ECdist_v5(PC_A1, PC_B1, PC_C1, PC_D1, normalization, odorNumber)
tmpd_origin=[0 0 0];





%% ----------------------------------
% AB
tmpVect1=PC_A1;
tmpVect2=PC_B1;

for ec=1:size(tmpVect1,1)
   tmpd_current1=tmpVect1(ec,:);
   tmpd_current2=tmpVect2(ec,:);
   X1X2(ec,1)=pdist([tmpd_current1;tmpd_current2]);
end

ECdist.AB1=X1X2;

if ~(odorNumber ==2)

% AC
tmpVect1=PC_A1;
tmpVect2=PC_C1;

for ec=1:size(tmpVect1,1)
   tmpd_current1=tmpVect1(ec,:);
   tmpd_current2=tmpVect2(ec,:);
   X1X2(ec,1)=pdist([tmpd_current1;tmpd_current2]);
end

ECdist.AC1=X1X2;


% AD
tmpVect1=PC_A1;
tmpVect2=PC_D1;

for ec=1:size(tmpVect1,1)
   tmpd_current1=tmpVect1(ec,:);
   tmpd_current2=tmpVect2(ec,:);
   X1X2(ec,1)=pdist([tmpd_current1;tmpd_current2]);
end

ECdist.AD1=X1X2;


% BC
tmpVect1=PC_B1;
tmpVect2=PC_C1;

for ec=1:size(tmpVect1,1)
   tmpd_current1=tmpVect1(ec,:);
   tmpd_current2=tmpVect2(ec,:);
   X1X2(ec,1)=pdist([tmpd_current1;tmpd_current2]);
end

ECdist.BC1=X1X2;


% BD
tmpVect1=PC_B1;
tmpVect2=PC_D1;

for ec=1:size(tmpVect1,1)
   tmpd_current1=tmpVect1(ec,:);
   tmpd_current2=tmpVect2(ec,:);
   X1X2(ec,1)=pdist([tmpd_current1;tmpd_current2]);
end

ECdist.BD1=X1X2;

% CD
tmpVect1=PC_C1;
tmpVect2=PC_D1;

for ec=1:size(tmpVect1,1)
   tmpd_current1=tmpVect1(ec,:);
   tmpd_current2=tmpVect2(ec,:);
   X1X2(ec,1)=pdist([tmpd_current1;tmpd_current2]);
end

ECdist.CD1=X1X2;



%% z-score Normalization by distances in Pre-stimulus period

if normalization ==1

preDist= [ECdist.AB1(1:20) ; ECdist.AC1(1:20)  ; ECdist.AD1(1:20)  ; ECdist.BC1(1:20)  ; ECdist.BD1(1:20)  ; ECdist.CD1(1:20)];

preDist_mean = mean(preDist);
preDist_std = std(preDist);


% ECdist.AB1=(ECdist.AB1 - preDist_mean)/preDist_std; 
% ECdist.AC1=(ECdist.AC1 - preDist_mean)/preDist_std;
% ECdist.AD1=(ECdist.AD1 - preDist_mean)/preDist_std;
% ECdist.BC1=(ECdist.BC1 - preDist_mean)/preDist_std;
% ECdist.BD1=(ECdist.BD1 - preDist_mean)/preDist_std;
% ECdist.CD1=(ECdist.CD1 - preDist_mean)/preDist_std;

ECdist.AB1=(ECdist.AB1 / preDist_mean); 
ECdist.AC1=(ECdist.AC1 / preDist_mean);
ECdist.AD1=(ECdist.AD1 / preDist_mean);
ECdist.BC1=(ECdist.BC1 / preDist_mean);
ECdist.BD1=(ECdist.BD1 / preDist_mean);
ECdist.CD1=(ECdist.CD1 / preDist_mean);


end

% ECdist.AB1=abs(ECdist.AB1); 
% ECdist.AC1=abs(ECdist.AC1);
% ECdist.AD1=abs(ECdist.AD1);
% ECdist.BC1=abs(ECdist.BC1);
% ECdist.BD1=abs(ECdist.BD1);
% ECdist.CD1=abs(ECdist.CD1);

end

end

