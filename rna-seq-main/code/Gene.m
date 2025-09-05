%% Genetic diversity analysis
%Check for unexpressed genes and delete them
n=all(HSMM==0,2)
HSMM(n,:)=[]
%Gene screening
m=size(HSMM,1)
for i=1:m
    e=find(HSMM(i,:)~=0)
    N(i)=length(e)
end
%Look for genes that are expressed in at least 250 cells
M=find(N>250)
Mean=zeros(m,1)
D=zeros(m,1)
CV=zeros(m,1)
m=size(M,2)
for i=1:m
    if any( i==M(:))
   Mean(i)=mean(HSMM(M(i),:))
   D(i)=std(HSMM(M(i),:))
   CV(i)=mean(i)/D(i)
    end
end
%Look for the top ten genes that are most different
m=size(M,2)
for i=1:m
    if CV(i)>10&&Mean(i)~=0
        Different(i)=i
    end
end
DIFFERENT=find(CV>10)

%% The expression of COL1A1 gene
%The original chronological gene expression was categorized by hours
gscatter(PCA1(:,1),PCA1(:,2),Type)
hold on
xlabel('CoL1A1');
title('CoL1A1 Gene Expression Change Map Arranged in Real Time');
%0/24/48/72h Line chart (four points represent the mean of each group)
plot(Psedotime(1,:),Psedotime(2,:))
%Sort pseudo time by hours
gscatter(Presedotime,COL1A1(:,1),Hours)
hold on
xlabel('CoL1A1');
title('CoL1A1 Gene Expression Change Map Arranged in Psedotime');
plot(TIME(1,:),TIME(2,:))
%According to Type classification pseudo time sort differentiation path separately marked
gscatter(Presedotime,COL1A1(:,1),Type)
hold on
xlabel('Psedotime');
title('CoL1A1 Gene Expression Change Map Arranged in Psedotime');
plot(A(1,:),A(2,:))%First differentiation pathway
plot(B(1,:),B(2,:))%Second differentiation pathway
legend('Type1','Type2','Type3','Type4','Type5','cell differentiation pathway 1','cell differentiation pathway 2')

%% The expression of ATG12 gene
%The original chronological gene expression was categorized by hours
gscatter(Gene,ATG12(:,1),Hours)
hold on
xlabel('ATG12');
title('ATG12 Gene Expression Change Map Arranged in Real Time');
%0/24/48/72h Line chart (four points represent the mean of each group)
plot(Psedotime(1,:),Psedotime(3,:))
%According to Type classification pseudo time sort differentiation path separately marked
gscatter(Presedotime,ATG12(:,1),Type)
hold on
xlabel('Psedotime');
title('ATG12 Gene Expression Change Map Arranged in Psedotime');
plot(A(1,:),A(2,:))
plot(B(1,:),B(2,:))
legend('Type1','Type2','Type3','Type4','Type5','cell differentiation pathway 1','cell differentiation pathway 2')

%% The expression of SURF4 gene
%The original chronological gene expression was categorized by hours
gscatter(Gene,SURF4(:,1),Hours)
hold on
xlabel('SURF4');
title('SURF4 Gene Expression Change Map Arranged in Real Time');
%0/24/48/72h Line chart (four points represent the mean of each group)
plot(Psedotime(1,:),Psedotime(4,:))
%According to Type classification pseudo time sort differentiation path separately marked
gscatter(Presedotime,SURF4(:,1),Type)
hold on
xlabel('Psedotime');
title('SURF4 Gene Expression Change Map Arranged in Psedotime');
plot(A(1,:),A(2,:))%First differentiation pathway
plot(B(1,:),B(2,:))%Second differentiation pathway
legend('Type1','Type2','Type3','Type4','Type5','cell differentiation pathway 1','cell differentiation pathway 2')

%%Heat map of genetic change
%The expression matrix of ten genes was normalized
Gene=mapminmax(Gene,0,1)
h = heatmap(Gene);
h=imagesc(Gene)
%Draw scatter plots and lines in time
gscatter(X,Gene1(:,4),Hours)
hold on
plot(TIME(1,:),TIME(2,:))
%Draw scatter plots and lines in pseudo-time
gscatter(X,Gene(:,4),Hours)
hold on
plot(Psedotime(1,:),Psedotime(2,:))