%% Draw the initial graph, the input is the two-dimensional coordinates after dimensionality reduction, 
%% a matrix with n rows and 2 columns, and the output is the initial graph
%Calculate the Mahalanobis distance between cell nodes
[m,n] = size(PCA1);
Dis = ones(m,m);
Cov = cov(PCA1);
for i=1:m
for j1=1:m
D(i,j1)=((PCA1(i,:)-PCA1(j1,:))*inv(Cov)*(PCA1(i,:)-PCA1(j1,:))')^0.5;
end
end
D
%Create a graph
D1=sparse(D)
%Find the smallest tree and get the adjacency matrix
MST=mst(D1)
%Convert a weighted adjacency matrix to an unweighted adjacency matrix
MST1=find(MST)
for i=MST1				
    MST(i)=1;		
end
%The generated raw diagram
gplot(full(MST),PCA1)
hold on


%% Calculate some basic properties of total edge number, 
%% summary nodes, clustering coefficient, betweenness centrality;
%Calculate the total number of edges in the graph
n = num_edges(MST)
%Calculate the total number of nodes in the diagram
m = num_vertices(MST);
%The betweenness centrality of the compute node and each edge, 
%where the centrality of the center storage node; E stores the centrality of each edge
%The number of times a node serves as the shortest bridge between the other two nodes. 
%The higher the number of times a node acts as a mediator, the greater its betweenness centrality.
[center,E] = betweenness_centrality(MST)

%% Find the degree and distribution curve of each node in the network diagram
%MST -- adjacency matrix after preliminary deletion of nodes
%DeD - The degree distribution of each node in a network graph
%aver_DeD -- Average degree of the network graph
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
s1=sum(MST(:,i))
s2=sum(MST(i,:))
   if s1>=s2
    DeD(i)=s1
  else
    DeD(i)=s2
   end
end
aver_DeD=mean(DeD);
if sum(DeD)==0
disp('The network diagram consists only of isolated points');
return;
else
figure;
bar([1:m],DeD);
xlabel('n');
ylabel('K');
title('Moderate size distrbution of network nodes');
end
figure;
M=max(DeD);
%The maximum degree of nodes in the network diagram is M, 
%but the existence of nodes whose degree is 0 should be considered at the same time
for i=1:M+1; 
N_DeD(i)=length(find(DeD==i-1));
end
P_DeD=zeros(1,M+1);
P_DeD(:)=N_DeD(:)./sum(N_DeD);
bar([0:M],P_DeD,'r');
xlabel('K');
ylabel('P(K)');
title('Probability distribution map of node degree');

%% Graph compression
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)

%% Preprocessing nodes with poor mediation
center= betweenness_centrality(MST)
a=find(center==0)
while a>=10
   i=a
   MST([i],:)=0
   MST(:,[i])=0
   PCA1([i],:)=0
   center= betweenness_centrality(MST)
   a=find(center==0)
end
%Update node degree
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
DeD(i)=sum(MST(i,:))
end
aver_DeD=mean(DeD);
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)
%Drawing STEP1
gplot(full(MST),PCA1)
hold on
title('Trajectory compression step1')
xlabel('PCA1')
ylabel('PCA2')
legend('Original trajectory','Post-preprocessing trajectary')

%% Vertices with a degree of 1
m=size(MST,2)
for i=1:m
    if any( i==DeD1(:))
        MST(i,:)=0
        MST(:,i)=0
    else
    end
end
gplot(full(MST),PCA1)
hold on
title('Trajectory compression step1.1')
xlabel('PCA1')
ylabel('PCA2')
legend('Original trajectory','Post-preprocessing trajectary')
%Update node degree
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
DeD(i)=sum(MST(i,:))
end
aver_DeD=mean(DeD);
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)


%% Vertices with a degree of 2
%% Calculate the Angle between the points
%% Calculate the degree of Angle between degree 2 and the adjacent point
A=MST
pos=zeros(1,m)
angle=zeros(1,m)
realangle=zeros(1,m)
for i=1:m
   if any( i == DeD2(:))
      j1=find(A(i,:)~=0)
      %定位 度为2的节点
      x1 =PCA1(j1(1,1),1);y1 =PCA1(j1(1,1),2);
      x2 =PCA1(i,1);y2 =PCA1(i,2);
      x3=PCA1(j1(1,2),1);y3 = PCA1(j1(1,2),2);
      a2 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
      b2 = (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2);
      c2 = (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3);
      a = sqrt(a2);
      b = sqrt(b2);
      c = sqrt(c2);
     pos(i) = (a2+b2-c2)/(2*a*b);    %Let's figure out the cosine
     angle(i) = acos(pos(i));         %Let's replace cosine with radians
     realangle(i) = angle(i)*180/pi;   %The radians are converted to the angles
   else
   end
end
%% Vertices with a degree of 2 are processed, and weight updates are performed
m=size(A,2)
for i=1:m
   if any( i == DeD2(:))&&realangle(i)<130
      j1=find(A(i,:)~=0)
      %The preprocessed path does not contain loops, 
      %so for nodes with degree 2, they will randomly merge to a nearby vertex
      %If there is no edge, add an edge, and if there is an edge, update it according to the rule weight
          if A(j1(1,1),j1(1,2))~=0
              A(j1(1,1),j1(1,2))=A(j1(1,1),j1(1,2))+0.5*A(i,j1(1,1))*A(i,j1(1,2))
          else
              A(j1(1,1),j1(1,2))=1
              A(j1(1,2),j1(1,1))=1
          end
          A(i,:)=0
          A(:,i)=0
          %And delete the original branch of node i
   else  
   end
end
%View iteration diagram
gplot(full(A),PCA1)
title('Trajectory compression step2')
xlabel('PCA1')
ylabel('PCA2')
legend('Original trajectory','Post-preprocessing trajectary')
hold on
%Update node degree
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
DeD(i)=sum(MST(i,:))
end
aver_DeD=mean(DeD);
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)

%% Step2.1
m=size(A)
for i=1:m
    if any( i==DeD1(:))
        A(i,:)=0
        A(:,i)=0
    else
    end
end
%Update node degree
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
DeD(i)=sum(MST(i,:))
end
aver_DeD=mean(DeD);
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)
%View iteration diagram
gplot(full(A),PCA1)
title('Trajectory compression step2.1')
xlabel('PCA1')
ylabel('PCA2')
legend('Original trajectory','Post-preprocessing trajectary')
hold on

%% Vertices with a degree of 3
% The influence of the movement of each node with degree 3 on Q is calculated by modular merging method
%Calculate ∑i, which is the sum of adjacent weights of each node. 
%The initial state is the case of degrees, because each branch of the undirected graph has a weight of 1
sum_i=DeD
%Calculate the change in Q value when any two nodes i merge into j
m3=size(DeD3,2)
N=sum((sum(A))')
Q_changevalue=zeros(m,m)
%j is the three nodes adjacent to the node of degree 3. 
%The change of Q value when the point i is moved to the point j is calculated respectively
for i=1:m
   if any( i == DeD3(:))
      j1=find(A(i,:)~=0)
      for c=1:3
          Inside=(sum_i(1,j1(1,c))+A(i,j1(1,c)))/(2*N)
          Outside=(sum_i(1,i)/(2*N))^2
          Q_changevalue(i,j1(1,c))=Inside-Outside%Q是一个上三角矩阵
      end
   else    
   end
end
Power=A
for i=1:m
    if any( i==DeD3(:))
        j1=find(Power(i,:)~=0)
        %If the modularity coefficient of the first 
        %node increases the most among the three adjacent nodes, combine node i and j(1,1).
        if Q_changevalue(i,j1(1,1))>= Q_changevalue(i,j1(1,2))&&Q_changevalue(i,j1(1,1))>= Q_changevalue(i,j1(1,3))
            Power(j1(1,1),j1(1,2))=Power(j1(1,1),j1(1,2))+Power(i,j1(1,1))
            Power(j1(1,1),j1(1,3))=Power(j1(1,1),j1(1,3))+Power(i,j1(1,3))
        else
        end
        %If the modularity coefficient of the second 
        %node increases the most among the three adjacent nodes, i node and j(1,2) are merged.
        if Q_changevalue(i,j1(1,2))> Q_changevalue(i,j1(1,1))&&Q_changevalue(i,j1(1,2))>= Q_changevalue(i,j1(1,3))
           Power(j1(1,2),j1(1,1))=Power(j1(1,2),j1(1,1))+Power(i,j1(1,1))
           Power(j1(1,2),j1(1,3))=Power(j1(1,2),j1(1,3))+Power(i,j1(1,3))
        else
        end
        %If the modularity coefficient of the third 
        %node increases the most among the three adjacent nodes, i node and j(1,3) are merged.
        if Q_changevalue(i,j1(1,3))> Q_changevalue(i,j1(1,1))&&Q_changevalue(i,j1(1,3))> Q_changevalue(i,j1(1,2))
           Power(j1(1,3),j1(1,1))=Power(j1(1,3),j1(1,1))+Power(i,j1(1,1))
           Power(j1(1,3),j1(1,2))=Power(j1(1,3),j1(1,2))+Power(i,j1(1,2))
        else
        end
    else
    end
end

for i=1:m
    if any( i==DeD3(:))
        Power(i,:)=0
        Power(:,i)=0
    else
    end
end
gplot(Power,PCA1)
hold on

%% Update Q and iterate with node case
%Update node degree
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
DeD(i)=sum(MST(i,:))
end
aver_DeD=mean(DeD);
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)

%The influence of the movement of each node with degree 3 on Q is calculated by modular merging method
%Calculate ∑i, which is the sum of adjacent weights of each node. 
%The initial state is the case of degrees, because each branch of the undirected graph has a weight of 1
sum_i=DeD
%Calculate the change in Q value when any two nodes i merge into j
m3=size(DeD3,2)
N=sum((sum(Power))')
Q_changevalue=zeros(m,m)
%j is the three nodes adjacent to the node of degree 3. 
%The change of Q value when the point i is moved to the point j is calculated respectively
for i=1:m
   if any( i == DeD3(:))
      j1=find(Power(i,:)~=0)
      for c=1:3
          Inside=(sum_i(1,j1(1,c))+Power(i,j1(1,c)))/(2*N)
          Outside=(sum_i(1,i)/(2*N))^2
          Q_changevalue(i,j1(1,c))=Inside-Outside
          %Q is an upper triangular matrix
      end
   else    
   end
end
title('Trajectory compression step3')
xlabel('PCA1')
ylabel('PCA2')
legend('Original trajectory','Post-preprocessing trajectary')
hold on
m=size(Power)
for i=1:m
    if any( i==DeD1(:))
        Power(i,:)=0
        Power(:,i)=0
    else
    end
end
gplot(Power,PCA1)
title('Trajectory compression step3.1')
xlabel('PCA1')
ylabel('PCA2')
legend('Original trajectory','Post-preprocessing trajectary')
hold on

%% Contraction step (if there are vertices of more than 4 degrees)
Zip=Power
sum_i=Power
%Calculate the change in Q value when any two nodes i merge into j
m3=size(DeD4,2)
N=sum((sum(Power))')
Q_changevalue=zeros(m,m)
%j is the three nodes adjacent to the node of degree 3. 
%The change of Q value when the point i is moved to the point j is calculated respectively
for i=1:m
   if any( i == DeD4(:))
      j1=find(Power(i,:)~=0)
      for c=1:4
          Inside=(sum_i(1,j1(1,c))+MST(i,j1(1,c)))/(2*N)
          Outside=(sum_i(1,i)/(2*N))^2
          Q_changevalue(i,j1(1,c))=Inside-Outside
          %Q is an upper triangular matrix
      end

   else    
   end
end
for i=1:m
    if any( i==DeD4(:))
        j1=find(Power(i,:)~=0)
        %If the modularity coefficient of the first 
        %node increases the most among the four adjacent nodes, i node and j(1,1) are merged.
        if Q_changevalue(i,j1(1,1))>= Q_changevalue(i,j1(1,2))&&Q_changevalue(i,j1(1,1))>= Q_changevalue(i,j1(1,3))&&Q_changevalue(i,j1(1,1))>= Q_changevalue(i,j1(1,4))
            Power(j1(1,1),j1(1,2))=Power(j1(1,1),j1(1,2))+Power(i,j1(1,1))
            Power(j1(1,1),j1(1,3))=Power(j1(1,1),j1(1,3))+Power(i,j1(1,3))
            Power(j1(1,1),j1(1,4))=Power(j1(1,1),j1(1,4))+Power(i,j1(1,4))
        else
        end
        %If the modularity coefficient of the second 
        %node increases the most among the four adjacent nodes, i node and j(1,2) are merged.
        if Q_changevalue(i,j1(1,2))> Q_changevalue(i,j1(1,1))&&Q_changevalue(i,j1(1,2))>= Q_changevalue(i,j1(1,3))&&Q_changevalue(i,j1(1,2))>= Q_changevalue(i,j1(1,4))
           Power(j1(1,2),j1(1,1))=Power(j1(1,2),j1(1,1))+Power(i,j1(1,1))
           Power(j1(1,2),j1(1,3))=Power(j1(1,2),j1(1,3))+Power(i,j1(1,3))
           Power(j1(1,2),j1(1,4))=Power(j1(1,2),j1(1,4))+Power(i,j1(1,4))
        else
        end
        %If the modularity coefficient of the third 
        %node increases the most among the four adjacent nodes, i node and j(1,3) are merged.
        if Q_changevalue(i,j1(1,3))> Q_changevalue(i,j1(1,1))&&Q_changevalue(i,j1(1,3))> Q_changevalue(i,j1(1,2))&&Q_changevalue(i,j1(1,3))>= Q_changevalue(i,j1(1,4))
           Power(j1(1,3),j1(1,1))=Power(j1(1,3),j1(1,1))+Power(i,j1(1,1))
           Power(j1(1,3),j1(1,2))=Power(j1(1,3),j1(1,2))+Power(i,j1(1,2))
           Power(j1(1,3),j1(1,4))=Power(j1(1,3),j1(1,4))+Power(i,j1(1,4))
        else
        end
        %If the modularity coefficient of the fourth 
        %node increases the most among the four adjacent nodes, i node and j(1,4) are merged.
        if Q_changevalue(i,j1(1,4))> Q_changevalue(i,j1(1,1))&&Q_changevalue(i,j1(1,4))> Q_changevalue(i,j1(1,2))&&Q_changevalue(i,j1(1,4))> Q_changevalue(i,j1(1,3))
           Power(j1(1,4),j1(1,1))=Power(j1(1,4),j1(1,1))+Power(i,j1(1,1))
           Power(j1(1,4),j1(1,2))=Power(j1(1,4),j1(1,2))+Power(i,j1(1,2))
           Power(j1(1,4),j1(1,3))=Power(j1(1,4),j1(1,3))+Power(i,j1(1,3))
        else
        end
    else
    end
end

%% Update Q and iterate with node case
%Update node degree
m=size(MST,2);
DeD=zeros(1,m);
for i=1:m
DeD(i)=sum(MST(i,:))
end
aver_DeD=mean(DeD);
%View the node whose degree is 1/2/3/4
DeD1=find(DeD==1)
DeD2=find(DeD==2)
DeD3=find(DeD==3)
DeD4=find(DeD==4)

%Getting rid of the vertices with a degree of 4
for i=1:m
    if any( i==DeD4(:))
        Power(i,:)=0
        Power(:,i)=0
    else
    end
end
gplot(Power,PCA1)
hold on
%A node of degree 4 is reduced to a node of degree 3 after processing

%% Draw the final graph compression result and find the key nodes
title('Pseudo-time trajectory of HSMM dataset based on graph compression results')
xlabel('PCA1')
ylabel('PCA2')
legend('Post-preprocessing trajectary','Original trajectory')
hold on

%% Calculate the pseudo time of each node in the compressed graph
%Calculate the distance between the remaining nodes, here using the Dijstra algorithm
%Calculate the dijkstra distance between the four important nodes and the starting point
%Important nodes are nodes whose degree is greater than 1 and the starting cell. Four nodes are used as examples here
Persedotime=dijkstra_sp(MST,58)%Starting cell
%Calculate the dijkstra distance from four important nodes to each point
PersedotimeA1=dijkstra_sp(MST,91)%Branch
PersedotimeA2=dijkstra_sp(MST,158)%branch
PersedotimeA3=dijkstra_sp(MST,264)%branch
PersedotimeA4=dijkstra_sp(MST,136)%branch
%All cell points were categorized to calculate which branch point was closer
m=size(MST)
for i=1:m
        if PersedotimeA1(i)<PersedotimeA2(i)&&PersedotimeA1(i)<PersedotimeA3(i)&&PersedotimeA1(i)<PersedotimeA4(i)&&PersedotimeA1(i)<Persedotime(i)
        Type(i)=1
    else
    end
        if PersedotimeA2(i)<PersedotimeA1(i)&&PersedotimeA2(i)<PersedotimeA3(i)&&PersedotimeA2(i)<PersedotimeA4(i)&&PersedotimeA2(i)<Persedotime(i)
        Type(i)=2
    else
    end
        if PersedotimeA3(i)<PersedotimeA1(i)&&PersedotimeA3(i)<PersedotimeA2(i)&&PersedotimeA3(i)<PersedotimeA4(i)&&PersedotimeA3(i)<Persedotime(i)
        Type(i)=3
    else
    end
        if PersedotimeA4(i)<PersedotimeA1(i)&&PersedotimeA4(i)<PersedotimeA2(i)&&PersedotimeA4(i)<PersedotimeA3(i)&&PersedotimeA4(i)<Persedotime(i)
        Type(i)=4
    else
    end
    if  Persedotime(i)<PersedotimeA1(i)&&Persedotime(i)<PersedotimeA2(i)&&Persedotime(i)<PersedotimeA3(i)&&Persedotime(i)<PersedotimeA4(i)
        Type(i)=0
    else
    end
end
%The pseudo-time of each cell point was calculated
for i=1:m
    if Type(i)==1
    PERSEDOTIME(i)=Persedotime(91)+PersedotimeA1(i)
    end
        if Type(i)==2
           PERSEDOTIME(i)=Persedotime(158)+PersedotimeA2(i)
        end
        if Type(i)==3
           PERSEDOTIME(i)=Persedotime(264)+PersedotimeA3(i)
        end
        if Type(i)==4
           PERSEDOTIME(i)=Persedotime(136)+PersedotimeA4(i)
        end
         if Type(i)==0
           PERSEDOTIME(i)=Persedotime(i)
    end
end
gscatter(PCA1(:,1),PCA1(:,2),Type)
hold on
gplot(MST,PCA1)
title('Cell differentiation diagram')
xlabel('PCA1')
ylabel('PCA2')
legend('0','1','2','3','4','MST')

%% Typing(for example:3)
for i=1:m
        if PersedotimeA1(i)<PersedotimeA2(i)
        Type(i)=1
        else
        end
        if PersedotimeA1(i)==PersedotimeA2(i)
        Type(i)=2
        else
        end
        if PersedotimeA1(i)>PersedotimeA2(i)
        Type(i)=3
        else
        end
end
