        clear        
        close all
        cd('D:\Project\Neuron EV\Data');%Set the main folder directory.
        mainfolder=cd;
        samplecd=fullfile(mainfolder, 'AD CSF');%Set the directory for EV lists of individual patients.
        controlcd=fullfile(mainfolder, 'HC CSF');%Set the directory for EV lists of individual controls.
        combinedcd=fullfile(mainfolder, 'CombinedFiles');%Set the directory for combined EV lists.
        cd(combinedcd);
        fstruct=dir('*.txt');
        A="Running. The files that will be analyzed are:";
        for cnt=1:length(fstruct)
          A=[A; fstruct(cnt).name];
        end
        fprintf( '%s\n', A);
        EVlist=[];
        for cnt=1:length(fstruct)
          fnam=fstruct(cnt).name;
          rawdata=importdata(fnam);
          rawdata(find(rawdata<0))=0;
          positivedata=rawdata(find(sum(rawdata,2)>0),:);
          Size(cnt) = length(positivedata);
          EVlist=[EVlist; positivedata];
        end
        for i=1:sum(Size)
        for m=1:length(fstruct)-1
            if i<Size(1)+1
                group(i)=1;%Generate IDs to correlate EVs with folder.
            elseif i>sum(Size(1:m)) && i<sum(Size(1:m+1))+1
                group(i)=m+1;
            end
        end
        end
        Nmarker=min(size(EVlist));
        EVstandard1=zeros(Nmarker,nchoosek(Nmarker,1)*3);%Defining possible EV subtypes.
        Combo=nchoosek(1:Nmarker,1);
        for i=1:nchoosek(Nmarker,1)
            for j=1:3
        EVstandard1(Combo(i,1),(i-1)*3+j)=j;
            end
        end
        if Nmarker>1
        EVstandard2=zeros(Nmarker,nchoosek(Nmarker,2)*9);
        Combo=nchoosek(1:Nmarker,2);
        for i=1:nchoosek(Nmarker,2)
            for j=1:3
                for k=1:3
        EVstandard2(Combo(i,:),(i-1)*9+(j-1)*3+k)=[j k];
                end
            end
        end
        else
        EVstandard2=[];
        end
        if Nmarker>2
        EVstandard3=zeros(Nmarker,nchoosek(Nmarker,3)*27);
        Combo=nchoosek(1:Nmarker,3);
        for i=1:nchoosek(Nmarker,3)
            for j=1:3
                for k=1:3
                    for l=1:3
        EVstandard3(Combo(i,:),(i-1)*27+(j-1)*9+(k-1)*3+l)=[j k l];
                    end
                end
            end
        end
        else
        EVstandard3=[]; 
        end
        if Nmarker>3
        EVstandard4=zeros(Nmarker,nchoosek(Nmarker,4)*81);
        Combo=nchoosek(1:Nmarker,4);
        for i=1:nchoosek(Nmarker,4)
            for j=1:3
                for k=1:3
                    for l=1:3
                        for m=1:3
        EVstandard4(Combo(i,:),(i-1)*81+(j-1)*27+(k-1)*9+(l-1)*3+m)=[j k l m];
                        end
                    end
                end
            end
        end
        else
        EVstandard4=[];   
        end
        EVstandard1=[EVstandard1 EVstandard2 EVstandard3 EVstandard4]; %No more than 4 marker colocalization in this example.
        EVstandardnorm=EVstandard1;
        low1=zeros(1,Nmarker);
        mid1=zeros(1,Nmarker);
        high1=zeros(1,Nmarker);
        low2=zeros(1,Nmarker);
        mid2=zeros(1,Nmarker);
        high2=zeros(1,Nmarker);
        for i=1:Nmarker
        sortI=sort(EVlist(:,i));
        postI=sortI(find(sortI));
        Ie(i,1)=0.5*postI(round(length(postI)/3));%Define centers of grids.
        Ie(i,2)=1.5*postI(round(length(postI)/3));
        Ie(i,3)=(2*postI(round(length(postI)*2/3))-1.5*postI(round(length(postI)/3)));
        EVstandard1(i,find(EVstandard1(i,:)==1))=Ie(i,1);
        EVstandard1(i,find(EVstandard1(i,:)==2))=Ie(i,2);
        EVstandard1(i,find(EVstandard1(i,:)==3))=Ie(i,3);
        low1(i)=0;
        low2(i)=postI(round(length(postI)/3));%Define grid boundaries.
        mid1(i)=postI(round(length(postI)/3));
        mid2(i)=postI(round(length(postI)*2/3));
        high1(i)=postI(round(length(postI)*2/3));
        high2(i)=postI(round(length(postI)*0.95));
        end
        EVstandard=EVstandard1';
        subcount=zeros(length(EVstandard),length(fstruct));%Establish grids based on EV subtypes.
        EVstandardunit=EVstandard;
        EVstandardunit(find(EVstandardunit))=1;
        EVlistunit=EVlist;
        EVlistunit(find(EVlistunit))=1;
        groupassign=zeros(length(EVlist),1);%Subtype assignment array for all EVs.
        parfor i=1:length(EVlist)
        EVprofile=EVlist(i,:);
        EVprofileunit=EVlistunit(i,:);
        EVsub=(EVstandardunit-EVprofileunit);
        index=find(any(EVsub')==0);
        if isempty(index)==0
        Eucdistance=sum((EVstandard(index,:)-EVprofile).^2,2);
        index2=find(Eucdistance==min(Eucdistance));
        groupassign(i)=min(index(index2));%Assign EV to subtypes.
        end
        end
        for i=1:length(EVstandard)
            for j=1:max(group)
                index= find(group==j);
            subcount(i,j)=length(find(groupassign(index)==i));%Count the population of each grid(subtype).
            end
        end
            cd(samplecd);
            fstruct=dir('*.txt');
            indexR=find((subcount(:,2)+subcount(:,1))./length(EVlist)>0 & abs((subcount(:,2))./length(find(group==2))-(subcount(:,1))./length(find(group==1)))>0);%Remove empty bins and bins with no change under disease condition.
            EVstandard=EVstandardnorm(:,indexR)';
            [m,n]=size(EVstandard);
           subgroup=zeros(length(fstruct),length(EVstandard));
           EVsl=length(EVstandard);
           parfor i=1:length(fstruct)
               fnam=fstruct(i).name;
               EVlisttemp=importdata(fnam);
               for j=1:EVsl
                   panel=EVstandard(j,:);
                   index1=find(panel==1);
                   index2=find(panel==2);
                   index3=find(panel==3);
                   index4=find(panel==0);
                   index5=find(panel>0);
                   panellow=zeros(1,Nmarker);
                   panelhigh=zeros(1,Nmarker);
                   panellow(1,index1)=low1(1,index1);
                   panelhigh(1,index1)=low2(1,index1);
                   panellow(1,index2)=mid1(1,index2);
                   panelhigh(1,index2)=mid2(1,index2);
                   panellow(1,index3)=high1(1,index3);
                   panelhigh(1,index3)=high2(1,index3);
                   m=0;
                   n=1;
                   while n<=length(EVlisttemp)
               if all(EVlisttemp(n,index4)==0) && all(EVlisttemp(n,index5)<panelhigh(1,index5)) && all(EVlisttemp(n,index5)>panellow(1,index5)) 
                   m=m+1;
               end
               n=n+1;
                   end
               subgroup(i,j)=m/length(EVlisttemp);%Count population of each subgroup in each sample.
               end
           end
           cd(controlcd);
           fstruct1=dir('*.txt');
           subgroup1=zeros(length(fstruct1),length(EVstandard));
           parfor i=1:length(fstruct1)
               fnam=fstruct1(i).name;
               EVlisttemp=importdata(fnam);
               for j=1:EVsl
                   panel=EVstandard(j,:);
                   index1=find(panel==1);
                   index2=find(panel==2);
                   index3=find(panel==3);
                   index4=find(panel==0);
                   index5=find(panel>0);
                   panellow=zeros(1,Nmarker);
                   panelhigh=zeros(1,Nmarker);
                   panellow(1,index1)=low1(1,index1);
                   panelhigh(1,index1)=low2(1,index1);
                   panellow(1,index2)=mid1(1,index2);
                   panelhigh(1,index2)=mid2(1,index2);
                   panellow(1,index3)=high1(1,index3);
                   panelhigh(1,index3)=high2(1,index3);
                   m=0;
                   n=1;
                   while n<=length(EVlisttemp)
               if all(EVlisttemp(n,index4)==0) && all(EVlisttemp(n,index5)<panelhigh(1,index5)) && all(EVlisttemp(n,index5)>panellow(1,index5)) 
                   m=m+1;
               end
               n=n+1;
                   end
               subgroup1(i,j)=m/length(EVlisttemp);%Count again for the control group.
               end
           end
          pvalues = mattest(subgroup', subgroup1');%Calculate P values.
          figure
          SigStructure=mavolcanoplot(subgroup',subgroup1',pvalues,'Plotonly',true,'LogTrans',true);%Generate volcano plot using P>=0.05, FC>2.
          EVIDs=SigStructure.GeneLabels;
          if length(EVIDs)<100 %In case the number of subgroup for optimization is very small, include more EVs.
          SigStructure=mavolcanoplot(subgroup',subgroup1',pvalues,'Plotonly',true,'LogTrans',true,'PCutoff',0.1);%Generate volcano plot using P>=0.1, FC>2..
          EVIDs=SigStructure.GeneLabels;
          end
          AUCorigin=zeros(1,length(EVIDs));
          AUCfinal=zeros(1,length(EVIDs));
          for i=1:length(EVIDs)
               panel=EVstandard(str2double(EVIDs(i)),:);
               index1=find(panel==1);
               index2=find(panel==2);
               index3=find(panel==3);
               uplimit=zeros(1,length(panel));
               lowlimit=zeros(1,length(panel));
               uplimit(1,index1)=low2(1,index1);
               uplimit(1,index2)=mid2(1,index2);
               uplimit(1,index3)=high2(1,index3);
               lowlimit(1,index1)=low1(1,index1);
               lowlimit(1,index2)=mid1(1,index2);
               lowlimit(1,index3)=high1(1,index3);
          end
          for i=1:length(EVIDs)
          fprintf('P-Values of EVID %4.2d: %4.2f\n',indexR(str2double(EVIDs(i))),SigStructure.PValues(i));%Print P-values.
          end
          for i=1:length(EVIDs)
          fprintf('Fold Changes of EVID %4.2d: %4.2f\n',indexR(str2double(EVIDs(i))),SigStructure.FoldChanges(i));%Print fold of change.
          end
          close figure 1
          parfor i=1:length(EVIDs) %Start a loop to optimize marker expression range.
               panel=EVstandard(str2double(EVIDs(i)),:);
               index1=find(panel==1);
               index2=find(panel==2);
               index3=find(panel==3);
               uplimit=zeros(1,length(panel));
               lowlimit=zeros(1,length(panel));
               uplimit(1,index1)=low2(1,index1);
               uplimit(1,index2)=mid2(1,index2);
               uplimit(1,index3)=high2(1,index3);
               lowlimit(1,index1)=low1(1,index1);
               lowlimit(1,index2)=mid1(1,index2);
               lowlimit(1,index3)=high1(1,index3);
               postindex=find(uplimit>0);
               zeroindex=find(uplimit==0);
               lowerthresh1=lowlimit(postindex);
               upperthresh1=uplimit(postindex);
               cd(samplecd);
               EVcount=zeros(length(fstruct),1);  
               for cnt=1:length(fstruct)
               fnam=fstruct(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               B=[];
               Comindx=[];
               while m <= length(EVlisttemp)
               B = find(EVlisttemp(m,postindex)>lowerthresh1 & EVlisttemp(m,postindex)<upperthresh1);
               if size(B,2) == length(postindex)
               Comindx(n+1)=m;
               n=n+1;              
               end
               m=m+1;
               end
               EVcount(cnt)=n/m*100;
               end
               cd(controlcd);
               EVcount1=zeros(length(fstruct1),1);  
               for cnt=1:length(fstruct1)
               fnam=fstruct1(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               B=[];
               Comindx=[];
               while m <= length(EVlisttemp)
               B = find(EVlisttemp(m,postindex)>lowerthresh1 & EVlisttemp(m,postindex)<upperthresh1);
               if size(B,2) == length(postindex)
               Comindx(n+1)=m;
               n=n+1;              
               end
               m=m+1;
               end
               EVcount1(cnt)=n/m*100;
               end
               EVcountsum=[EVcount1; EVcount];
               resp = (1:length(EVcountsum))'>length(EVcount1);
               [X,Y,T,AUCraw] = perfcurve(resp,EVcountsum,1);%Calculate unoptimized AUCs.
               if AUCraw<0.5
                   AUCraw=1-AUCraw;
               end
               if  round(AUCraw,2)>0.65
               AUCopt=AUCraw;
               lowlimit1=lowlimit;
               uplimit1=uplimit;
               AUCorigin(i)=AUCraw;
               end
               if round(AUCraw,2)>0.65
               if length(postindex)==1 %Optimize AUCs for EVs with single marker.
               for indx1=1:7
                   for indx2=1:7
               lowerthresh1=((indx1-4)*0.2+1)*lowlimit(postindex);
               upperthresh1=((indx2-4)*0.2+1)*uplimit(postindex);
               cd(samplecd);
               for cnt=1:length(fstruct)
               fnam=fstruct(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               B=[];
               while m <= length(EVlisttemp)
               B = find(EVlisttemp(m,postindex)>lowerthresh1 & EVlisttemp(m,postindex)<upperthresh1);
               if size(B,2) == length(postindex) 
               n=n+1;
               end
               m=m+1;
               end
               EVcount(cnt)=n/m*100;
               end
               cd(controlcd); 
               EVcount1=zeros(length(fstruct1),1);
               for cnt=1:length(fstruct1)
               fnam=fstruct1(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m <= length(EVlisttemp)
               B = find(EVlisttemp(m,postindex)>lowerthresh1 & EVlisttemp(m,postindex)<upperthresh1);
               if size(B,2) == length(postindex) 
               n=n+1;
               end
               m=m+1;
               end
               EVcount1(cnt)=n/m*100;
               end
               EVcountsum=[EVcount1;EVcount];
               resp = (1:length(EVcountsum))'>length(EVcount1);
               [X,Y,T,AUCtemp] = perfcurve(resp,EVcountsum,1);
               if AUCtemp<0.5
                   AUCtemp=1-AUCtemp;
               end
               if AUCtemp>AUCopt
               AUCopt=AUCtemp;
               lowlimit1(1,postindex)=lowerthresh1;
               uplimit1(1,postindex)=upperthresh1;
               end
                   end
               end
               end
               if length(postindex)==2 %Optimize AUCs for EVs with 2 markers.
               for indx1=1:7
                   for indx2=1:7
                       for indx3=1:7
                           for indx4=1:7
               lowerthresh1=((indx1-5)*0.2+1)*lowlimit(postindex(1));
               upperthresh1=((indx2-5)*0.2+1)*uplimit(postindex(1));
               lowerthresh2=((indx3-5)*0.2+1)*lowlimit(postindex(2));
               upperthresh2=((indx4-5)*0.2+1)*uplimit(postindex(2));
               cd(samplecd);
               for cnt=1:length(fstruct)
               fnam=fstruct(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m < length(EVlisttemp)
               if EVlisttemp(m,postindex(1))>lowerthresh1 && EVlisttemp(m,postindex(1))<upperthresh1 && EVlisttemp(m,postindex(2))>lowerthresh2 && EVlisttemp(m,postindex(2))<upperthresh2
               n=n+1;
               end              
               m=m+1;
               end
               EVcount(cnt)=n/m*100;
               end
               cd(controlcd);
               EVcount1=zeros(length(fstruct1),1);
               for cnt=1:length(fstruct1)
               fnam=fstruct1(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m <= length(EVlisttemp)
               if EVlisttemp(m,postindex(1))>lowerthresh1 && EVlisttemp(m,postindex(1))<upperthresh1 && EVlisttemp(m,postindex(2))>lowerthresh2 && EVlisttemp(m,postindex(2))<upperthresh2
               n=n+1;
               end              
               m=m+1;
               end
               EVcount1(cnt)=n/m*100;
               end
               EVcountsum=[EVcount1;EVcount];
               resp = (1:length(EVcountsum))'>length(EVcount1);
               [X,Y,T,AUCtemp] = perfcurve(resp,EVcountsum,1);
               if AUCtemp<0.5
                   AUCtemp=1-AUCtemp;
               end
               if AUCtemp>AUCopt
               AUCopt=AUCtemp;
               lowlimit1(1,postindex(1))=lowerthresh1;
               uplimit1(1,postindex(1))=upperthresh1;
               lowlimit1(1,postindex(2))=lowerthresh2;
               uplimit1(1,postindex(2))=upperthresh2;              
               end
                   end
               end
                   end
               end
               end
               if length(postindex)==3 %Optimize AUCs for EVs with 3 markers.
               for indx1=1:5
                   for indx2=1:5
                       for indx3=1:5
                           for indx4=1:5
                               for indx5=1:5
                                   for indx6=1:5
               lowerthresh1=((indx1-3)*0.25+1)*lowlimit(postindex(1));
               upperthresh1=((indx2-3)*0.25+1)*uplimit(postindex(1));
               lowerthresh2=((indx3-3)*0.25+1)*lowlimit(postindex(2));
               upperthresh2=((indx4-3)*0.25+1)*uplimit(postindex(2));
               lowerthresh3=((indx5-3)*0.25+1)*lowlimit(postindex(3));
               upperthresh3=((indx6-3)*0.25+1)*uplimit(postindex(3));
               cd(samplecd);
               for cnt=1:length(fstruct)
               fnam=fstruct(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m <= length(EVlisttemp)
               if EVlisttemp(m,postindex(1))>lowerthresh1 && EVlisttemp(m,postindex(1))<upperthresh1 && EVlisttemp(m,postindex(2))>lowerthresh2 && EVlisttemp(m,postindex(2))<upperthresh2 && EVlisttemp(m,postindex(3))>lowerthresh3 && EVlisttemp(m,postindex(3))<upperthresh3
               n=n+1;
               end              
               m=m+1;
               end
               EVcount(cnt)=n/m*100;
               end
               cd(controlcd); 
               EVcount1=zeros(length(fstruct1),1);
               for cnt=1:length(fstruct1)
               fnam=fstruct1(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m <= length(EVlisttemp)
               if EVlisttemp(m,postindex(1))>lowerthresh1 && EVlisttemp(m,postindex(1))<upperthresh1 && EVlisttemp(m,postindex(2))>lowerthresh2 && EVlisttemp(m,postindex(2))<upperthresh2 && EVlisttemp(m,postindex(3))>lowerthresh3 && EVlisttemp(m,postindex(3))<upperthresh3
               n=n+1;
               end              
               m=m+1;
               end
               EVcount1(cnt)=n/m*100;
               end
               EVcountsum=[EVcount1;EVcount];
               resp = (1:length(EVcountsum))'>length(EVcount1);
               [X,Y,T,AUCtemp] = perfcurve(resp,EVcountsum,1);
               if AUCtemp<0.5
                   AUCtemp=1-AUCtemp;
               end
               if AUCtemp>AUCopt
               AUCopt=AUCtemp;
               lowlimit1(1,postindex(1))=lowerthresh1;
               uplimit1(1,postindex(1))=upperthresh1;
               lowlimit1(1,postindex(2))=lowerthresh2;
               uplimit1(1,postindex(2))=upperthresh2;       
               lowlimit1(1,postindex(3))=lowerthresh3;
               uplimit1(1,postindex(3))=upperthresh3;  
               end
                   end
               end
                   end
               end
                   end
               end
               end
               if length(postindex)==4  %Optimize AUCs for EVs with 4 markers.
               for indx1=1:3
                   for indx2=1:3
                       for indx3=1:3
                           for indx4=1:3
                               for indx5=1:3
                                   for indx6=1:3
                                       for indx7=1:3
                                           for indx8=1:3
               lowerthresh1=((indx1-2)*0.5+1)*lowlimit(postindex(1));
               upperthresh1=((indx2-2)*0.5+1)*uplimit(postindex(1));
               lowerthresh2=((indx3-2)*0.5+1)*lowlimit(postindex(2));
               upperthresh2=((indx4-2)*0.5+1)*uplimit(postindex(2));
               lowerthresh3=((indx5-2)*0.5+1)*lowlimit(postindex(3));
               upperthresh3=((indx6-2)*0.5+1)*uplimit(postindex(3));
               lowerthresh4=((indx7-2)*0.5+1)*lowlimit(postindex(4));
               upperthresh4=((indx8-2)*0.5+1)*uplimit(postindex(4));
               cd(samplecd);
               for cnt=1:length(fstruct)
               fnam=fstruct(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m <= length(EVlisttemp)
               if EVlisttemp(m,postindex(1))>lowerthresh1 && EVlisttemp(m,postindex(1))<upperthresh1 && EVlisttemp(m,postindex(2))>lowerthresh2 && EVlisttemp(m,postindex(2))<upperthresh2 && EVlisttemp(m,postindex(3))>lowerthresh3 && EVlisttemp(m,postindex(3))<upperthresh3 && EVlisttemp(m,postindex(4))>lowerthresh4 && EVlisttemp(m,postindex(4))<upperthresh4
               n=n+1;
               end              
               m=m+1;
               end
               EVcount(cnt)=n/m*100;
               end
               cd(controlcd); 
               EVcount1=zeros(length(fstruct1),1);
               for cnt=1:length(fstruct1)
               fnam=fstruct1(cnt).name;
               EVlisttemp=importdata(fnam);
               n = 0;
               m = 1;
               while m <= length(EVlisttemp)
               if EVlisttemp(m,postindex(1))>lowerthresh1 && EVlisttemp(m,postindex(1))<upperthresh1 && EVlisttemp(m,postindex(2))>lowerthresh2 && EVlisttemp(m,postindex(2))<upperthresh2 && EVlisttemp(m,postindex(3))>lowerthresh3 && EVlisttemp(m,postindex(3))<upperthresh3 && EVlisttemp(m,postindex(4))>lowerthresh4 && EVlisttemp(m,postindex(4))<upperthresh4
               n=n+1;              
               end
               m=m+1;
               end
               EVcount1(cnt)=n/m*100;
               end
               EVcountsum=[EVcount1;EVcount];
               resp = (1:length(EVcountsum))'>length(EVcount1);
               [X,Y,T,AUCtemp] = perfcurve(resp,EVcountsum,1);
               if AUCtemp<0.5
                   AUCtemp=1-AUCtemp;
               end
               if AUCtemp>AUCopt
               AUCopt=AUCtemp;
               lowlimit1(1,postindex(1))=lowerthresh1;
               uplimit1(1,postindex(1))=upperthresh1;
               lowlimit1(1,postindex(2))=lowerthresh2;
               uplimit1(1,postindex(2))=upperthresh2;       
               lowlimit1(1,postindex(3))=lowerthresh3;
               uplimit1(1,postindex(3))=upperthresh3;
               lowlimit1(1,postindex(4))=lowerthresh4;
               uplimit1(1,postindex(4))=upperthresh4;  
               end
                   end
               end
                   end
               end
                   end
               end
                   end
               end
               end
               if round(AUCopt,2)>0.65
               fprintf('Expression Profile of EVID %4.2d: [%s]',indexR(str2double(EVIDs(i))),join(string(lowlimit1), ','));%Print the final expression profile.
               fprintf('-[%s]\n', join(string(uplimit1), ','));
               fprintf('EVID %4.2d, AUC range optimized: %4.2f\n',indexR(str2double(EVIDs(i))), AUCopt);
               fprintf('---------------------------------------------------------------\n');
               end
               AUCfinal(i)=AUCopt;
               end
          end