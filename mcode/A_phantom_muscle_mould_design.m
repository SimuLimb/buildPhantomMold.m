clear; close all; clc;

%% Build the phantom mould for muscle surface

%% Control parameters 
% Path names
projectFolder = fileparts(fileparts(mfilename('fullpath')));
projectFolder_sim_Amputation=fileparts('/Users/s2986149/Desktop/');
loadFolder=fullfile(projectFolder_sim_Amputation,'simulateAmputation.m','data','BodyParts3D','post'); 
saveFolder_stl=fullfile(projectFolder,'data','mould_stl'); 
saveFolder_mat=fullfile(projectFolder,'data','mould_mat'); 
saveNameGeom_tf='BodyParts3D_right_leg_transfemoral_amp_muscle_mould';

saveOn=1;
green = 1/255*[0, 100, 0];

%Geometric parameters
pointSpacing=5;
scaleFactor=[1.1 1.1 1.01];
scaleFactor1=[1.5 1.5 1];
CutHeight=30;
StandHeight=50;
select_amputation_case='tf';
switch select_amputation_case
     case 'tf'
        fileName ='BodyParts3D_right_leg_transfemoral_amp';
     case 'tt'
        fileName='BodyParts3D_right_leg_transtibial_amp';
end

%Load selected case of the processed amputation
fileName_mat=fullfile(loadFolder,[fileName,'.mat']);
model=load(fileName_mat);
F=model.FT_amp; %Faces
V=model.VT_amp; %Vertices
C=model.CT_amp; %Color label

switch select_amputation_case
    case 'tf'
        
        %Muscles
        F_muscle=F{3};
        V_muscle=V{3};
        [C,I]=max(V_muscle(:,3));

        %Skin
        F_skin=F{2};
        V_skin=V{2};

        %Femur
        F_femur=F{1};
        V_femur=V{1};

        %Use the most higher edge of the amputated muscle/smin+fat layer as the 
        %selected point for the slicing of the femur, to design the supported 
        %cylindrical bars going throuhg the femur and extended all the way to 
        %muscle and skin+fat layers laying on them.

        snapTolerance=mean(patchEdgeLengths(F_femur,V_femur))/100;
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P=V_muscle(I,:);%Point on plane corresponds to the higher edge of the muscle layer
        C=[];
        %Slicing surface 
        [F_sliced,V_sliced,~,logicSide,Eb]=triSurfSlice(F_femur,V_femur,C,P,n,snapTolerance);
        indCurveNow=edgeListToCurve(Eb);
        Vc=V_sliced(indCurveNow,:);

        %% Create and position cylinders geometry
        pointSpacingBeams=pointSpacing/2;
        inputStruct.cylRadius=5;
        inputStruct.numRadial=round((2*pi*inputStruct.cylRadius)./pointSpacingBeams);
        inputStruct.cylHeight=max(V_skin(:,1))-min(V_skin(:,1))+30;
        nh=round(inputStruct.cylHeight./pointSpacingBeams);
        nh=nh+double(iseven(nh));
        inputStruct.numHeight=nh;
        inputStruct.meshType='tri';
        inputStruct.closeOpt=1;
        % Derive patch data for a cylinder
        [F_cylinder,V_cylinder,C_cylinder]=patchcylinder(inputStruct);
        R1=euler2DCM([0.5*pi 0 0]);
        R2=euler2DCM([0 0.5*pi 0]);

        V_cylinder1=V_cylinder*R1+mean(Vc,1);
        V_cylinder1(:,2)=V_cylinder1(:,2)+8;
        V_cylinder2=V_cylinder*R2+mean(Vc,1);
        V_cylinder2(:,1)=V_cylinder1(:,2)-5;
        
        F_cylinder2=fliplr(F_cylinder);
        F_cylinder1=F_cylinder;
        
        % Visualize surface and landmarks
        cFigure; hold on;
        gpatch(F_femur,V_femur,'b','none',0.5);
        gpatch(F_muscle,V_muscle,'r','none',0.5);
        gpatch(F_skin,V_skin,green,'none',0.5);
        plotV(Vc,'k-','LineWidth',4);
       
        gpatch(F_cylinder1,V_cylinder1,'w','none',1);
        patchNormPlot(F_cylinder1,V_cylinder1);
        
        gpatch(F_cylinder2,V_cylinder2,'w','none',1);
        %patchNormPlot(F_cylinder2,V_cylinder2);
        
        axisGeom; axis off; camlight headlight;
        gdrawnow;

        % Find the intersection between cylinders and femur
        % Delete intersected points from the femur mesh

        %% First cylinder
        pointSpacing=mean(patchEdgeLengths(F_femur,V_femur));
        voxelSize=pointSpacing;
        [regionLabel_implant]=simplexImIntersect(F_femur,V_femur,[],V_cylinder1,voxelSize);
        logicOut_implant=~(~isnan(regionLabel_implant));
        logicFacesOut_implant=all(logicOut_implant(F_cylinder1),2);
        
        [regionLabel_bone]=simplexImIntersect(F_cylinder1,V_cylinder1,[],V_femur,voxelSize);
        logicOut_bone=isnan(regionLabel_bone);
        logicFacesOut_bone=all(logicOut_bone(F_femur),2);
        
        logicFacesOut_bone=triSurfLogicSharpFix(F_femur,logicFacesOut_bone,3);
        F1=F_femur(logicFacesOut_bone,:);
        V1=V_femur;
        [F1,V1]=patchCleanUnused(F1,V1);
        
        Eb=patchBoundary(F1,V1);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct); 
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve1=edgeListToCurve(E1);
        indCurve1=indCurve1(1:end-1);
        
        %[~,indKeep]=min(groupSize);
        logicKeep=G~=indKeep;
        E1=Eb(logicKeep,:);
        indCurve2=edgeListToCurve(E1);
        indCurve2=indCurve2(1:end-1);

        logicFacesIn_implant=all(~logicOut_implant(F_cylinder1),2);
        logicFacesIn_implant=triSurfLogicSharpFix(F_cylinder1,logicFacesIn_implant,3);
        F2=F_cylinder1(logicFacesIn_implant,:);
        V2=V_cylinder1;
        [F2,V2]=patchCleanUnused(F2,V2);
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [V2]=patchSmooth(F2,V2,[],cParSmooth);
        Eb=patchBoundary(F2,V2);
        
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E2=Eb(logicKeep,:);
        indCurve3=edgeListToCurve(E2);
        indCurve3=indCurve3(1:end-1);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E2=Eb(logicKeep,:);
        indCurve4=edgeListToCurve(E2);
        indCurve4=indCurve4(1:end-1);
        
        %Visualization
        cFigure; hold on;
        gpatch(F1,V1,'w','k',0.25);
        
        plotV(V1(indCurve2,:),'k.-','LineWidth',3);
        plotV(V2(indCurve3,:),'r.-','LineWidth',3);
        
        plotV(V1(indCurve1,:),'k.-','LineWidth',3);
        plotV(V2(indCurve4,:),'r.-','LineWidth',3);
        gpatch(F2,V2,'r','k',1)
        gpatch(F_cylinder1,V_cylinder1,'w','none',0.25)
        
%         
        axisGeom;axis off;
        camlight headlight;
        drawnow;

        %Find closer curve for lofting
        D1=min(minDist(V1(indCurve1,:),V2(indCurve3,:)));
        D2=min(minDist(V1(indCurve1,:),V2(indCurve4,:)));

        %Reorder curves
        if D2<D1
            [~,indMin]=minDist(V1(indCurve1(1),:),V2(indCurve4,:));
        
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            [Fh,Vh]=regionTriMesh3D({V1(indCurve1,:),V2(indCurve4,:)},pointSpacing,0,'natural');
        else
            [~,indMin]=minDist(V1(indCurve1(1),:),V2(indCurve3,:));
        
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            [Fh,Vh]=regionTriMesh3D({V1(indCurve1,:),V2(indCurve3,:)},pointSpacing,0,'natural');
        end
        
        %Find closer curve for lofting
        D1=min(minDist(V1(indCurve2,:),V2(indCurve3,:)));
        D2=min(minDist(V1(indCurve2,:),V2(indCurve4,:)));

        %Reorder curves
        if D2<D1
            [~,indMin]=minDist(V1(indCurve2(1),:),V2(indCurve4,:));            
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            [Fh1,Vh1]=regionTriMesh3D({V1(indCurve2,:),V2(indCurve4,:)},pointSpacing,0,'natural');
        else
           [~,indMin]=minDist(V1(indCurve2(1),:),V2(indCurve3,:));            
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            [Fh1,Vh1]=regionTriMesh3D({V1(indCurve2,:),V2(indCurve3,:)},pointSpacing,0,'natural');
        end
        
% 
%         %Visualization
%         cFigure; hold on;
%         gpatch(F1,V1,'w','k',1);
%         
%         plotV(V1(indCurve2,:),'g.-','LineWidth',5);
%         plotV(V2(indCurve3,:),'m.-','LineWidth',5);
%         
%         plotV(V1(indCurve1,:),'r.-','LineWidth',5);
%         plotV(V2(indCurve4,:),'c.-','LineWidth',5);
%         gpatch(F2,V2,'w','k',1)
% %         
%         axisGeom;
%         camlight headlight;
%         icolorbar;
%         drawnow;
        

        %Visualization
        cFigure; hold on;
        gpatch(F1,V1,'w','k',0.25);
        %patchNormPlot(F1,V1);

        gpatch(Fh,Vh,'g','k',1);
        %patchNormPlot(fliplr(Fh),Vh);

        gpatch(Fh1,Vh1,'r','k',1);
        %patchNormPlot(Fh1,Vh1);  

        axisGeom;axis off;
        camlight headlight;
        drawnow;

        %
        [Ff1,Vf1]=joinElementSets({F1,fliplr(Fh),Fh1},{V1,Vh,Vh1});
        [Ff1,Vf1]=patchCleanUnused(Ff1,Vf1);
        [Ff1,Vf1]=mergeVertices(Ff1,Vf1);

        
        Eb=patchBoundary(Ff1,Vf1);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve1=edgeListToCurve(E1);
        indCurve1=indCurve1(1:end-1);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve2=edgeListToCurve(E1);
        indCurve2=indCurve2(1:end-1);
        
        cFigure; hold on;
        gpatch(Ff1,Vf1,'w','k',1);
        patchNormPlot(Ff1,Vf1);  

        plotV(Vf1(indCurve1,:),'k.-','LineWidth',3);
        plotV(Vf1(indCurve2,:),'k.-','LineWidth',3);
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;

        %
        logicFacesOut_implant=all(logicOut_implant(F_cylinder1),2);
        logicFacesOut_implant=triSurfLogicSharpFix(F_cylinder1,logicFacesOut_implant,3);
        F1=F_cylinder1(logicFacesOut_implant,:);
        V1=V_cylinder1;
        [F1,V1]=patchCleanUnused(F1,V1);
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [V1]=patchSmooth(F1,V1,[],cParSmooth);
        Eb=patchBoundary(F1,V1);
        
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve3=edgeListToCurve(E1);
        indCurve3=indCurve3(1:end-1);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve4=edgeListToCurve(E1);
        indCurve4=indCurve4(1:end-1);
        
        cFigure; hold on;
        gpatch(Ff1,Vf1,'w','k',0.25);
        
        plotV(Vf1(indCurve1,:),'k.-','LineWidth',3);
        plotV(Vf1(indCurve2,:),'k.-','LineWidth',3);
        
        gpatch(F1,V1,'r','k',0.25)
        
        plotV(V1(indCurve3,:),'r.-','LineWidth',3);
        plotV(V1(indCurve4,:),'r.-','LineWidth',3);
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;
   
        D1=min(minDist(Vf1(indCurve1,:),V1(indCurve3,:)));
        D2=min(minDist(Vf1(indCurve1,:),V1(indCurve4,:)));
        
        
        if D2<D1
            [~,indMin]=minDist(Vf1(indCurve1(1),:),V1(indCurve4,:));            
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(indCurve4,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(flip(indCurve4),:)).^2,2)));
            if D2a<D1a
                indCurve4=flip(indCurve4);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh1,Vh1]=polyLoftLinear(Vf1(indCurve1,:),V1(indCurve4,:),cPar);
        else
           [~,indMin]=minDist(Vf1(indCurve1(1),:),V1(indCurve3,:));            
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(indCurve3,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(flip(indCurve3),:)).^2,2)));
            if D2a<D1a
                indCurve3=flip(indCurve3);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh1,Vh1]=polyLoftLinear(Vf1(indCurve1,:),V1(indCurve3,:),cPar);
        end

        D1=min(minDist(Vf1(indCurve2,:),V1(indCurve3,:)));
        D2=min(minDist(Vf1(indCurve2,:),V1(indCurve4,:)));

        if D2<D1
            [~,indMin]=minDist(Vf1(indCurve2(1),:),V1(indCurve4,:));            
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(indCurve4,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(flip(indCurve4),:)).^2,2)));
            if D2a<D1a
                indCurve4=flip(indCurve4);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh2,Vh2]=polyLoftLinear(Vf1(indCurve2,:),V1(indCurve4,:),cPar);
        else
           [~,indMin]=minDist(Vf1(indCurve2(1),:),V1(indCurve3,:));            
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(indCurve3,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(flip(indCurve3),:)).^2,2)));
            if D2a<D1a
                indCurve3=flip(indCurve3);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh2,Vh2]=polyLoftLinear(Vf1(indCurve2,:),V1(indCurve3,:),cPar);
        end
        
        
        cFigure; hold on;
        gpatch(Ff1,Vf1,'w','k',0.25);
        
        %plotV(Vf1(indCurve5,:),'g.-','LineWidth',5);
        %plotV(Vf1(indCurve6,:),'m.-','LineWidth',5);
        
        %plotV(V3(indCurve7,:),'c.-','LineWidth',5);
        %plotV(V3(indCurve8,:),'k.-','LineWidth',5);
        
        gpatch(fliplr(Fh1),Vh1,'r','k',1)
        %patchNormPlot(fliplr(Fh2),Vh2);  
        
        gpatch(fliplr(Fh2),Vh2,'g','k',1)
        %patchNormPlot(fliplr(Fh3),Vh3);  
        
        gpatch(F1,V1,'w','k',0.25)
        %patchNormPlot(F3,V3);
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;
  
        [F_bone,V_bone]=joinElementSets({Ff1,fliplr(Fh1),fliplr(Fh2),F1},{Vf1,Vh1,Vh2,V1});
        [F_bone,V_bone]=patchCleanUnused(F_bone,V_bone);
        [F_bone,V_bone]=mergeVertices(F_bone,V_bone);
        
        cFigure; hold on;
        gpatch(F_bone,V_bone,'w','k',1);
        patchNormPlot(F_bone,V_bone);  
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;

        %% Second cylinder
        pointSpacing=mean(patchEdgeLengths(F_bone,V_bone));
        voxelSize=pointSpacing;
        [regionLabel_implant]=simplexImIntersect(F_bone,V_bone,[],V_cylinder2,voxelSize);
        logicOut_implant=~(~isnan(regionLabel_implant));
        logicFacesOut_implant=all(logicOut_implant(F_cylinder2),2);
        
        [regionLabel_bone]=simplexImIntersect(F_cylinder2,V_cylinder2,[],V_bone,voxelSize);
        logicOut_bone=isnan(regionLabel_bone);
        logicFacesOut_bone=all(logicOut_bone(F_bone),2);
        
        logicFacesOut_bone=triSurfLogicSharpFix(F_bone,logicFacesOut_bone,3);
        F1=F_bone(logicFacesOut_bone,:);
        V1=V_bone;
        [F1,V1]=patchCleanUnused(F1,V1);
        
        Eb=patchBoundary(F1,V1);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct); 
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve1=edgeListToCurve(E1);
        indCurve1=indCurve1(1:end-1);
        
        %[~,indKeep]=min(groupSize);
        logicKeep=G~=indKeep;
        E1=Eb(logicKeep,:);
        indCurve2=edgeListToCurve(E1);
        indCurve2=indCurve2(1:end-1);

        logicFacesIn_implant=all(~logicOut_implant(F_cylinder2),2);
        logicFacesIn_implant=triSurfLogicSharpFix(F_cylinder2,logicFacesIn_implant,3);
        F2=F_cylinder2(logicFacesIn_implant,:);
        V2=V_cylinder2;
        [F2,V2]=patchCleanUnused(F2,V2);
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [V2]=patchSmooth(F2,V2,[],cParSmooth);
        Eb=patchBoundary(F2,V2);
        
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve3=edgeListToCurve(E1);
        indCurve3=indCurve3(1:end-1);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve4=edgeListToCurve(E1);
        indCurve4=indCurve4(1:end-1);
        
        %Find closer curve for lofting
        D1=min(minDist(V1(indCurve1,:),V2(indCurve3,:)));
        D2=min(minDist(V1(indCurve1,:),V2(indCurve4,:)));
        
        %Reorder curves
        if D2<D1
            [~,indMin]=minDist(V1(indCurve1(1),:),V2(indCurve4,:));
        
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            [Fh1,Vh1]=regionTriMesh3D({V1(indCurve1,:),V2(indCurve4,:)},pointSpacing,0,'natural');
        else
            [~,indMin]=minDist(V1(indCurve1(1),:),V2(indCurve3,:));
        
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            [Fh1,Vh1]=regionTriMesh3D({V1(indCurve1,:),V2(indCurve3,:)},pointSpacing,0,'natural');
        end
        
        D1=min(minDist(V1(indCurve2,:),V2(indCurve3,:)));
        D2=min(minDist(V1(indCurve2,:),V2(indCurve4,:)));

        %Reorder curves
        if D2<D1
            [~,indMin]=minDist(V1(indCurve2(1),:),V2(indCurve4,:));
        
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            [Fh2,Vh2]=regionTriMesh3D({V1(indCurve2,:),V2(indCurve4,:)},pointSpacing,0,'natural');
        else
            [~,indMin]=minDist(V1(indCurve2(1),:),V2(indCurve3,:));
        
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            [Fh2,Vh2]=regionTriMesh3D({V1(indCurve2,:),V2(indCurve3,:)},pointSpacing,0,'natural');
        end
        
        %Visualization
        cFigure; hold on;

         plotV(V1(indCurve1,:),'k.-','LineWidth',4);
         plotV(V1(indCurve2,:),'k.-','LineWidth',4);
%         
         plotV(V2(indCurve3,:),'r.-','LineWidth',4);
         plotV(V2(indCurve4,:),'r.-','LineWidth',4);
         
         gpatch(F1,V1,'w','k',0.25)
         gpatch(F2,V2,'r','k',0.25);
         
         gpatch(Fh1,Vh1,'g','k',1)
         patchNormPlot(Fh1,Vh1);
         
         gpatch(fliplr(Fh2),Vh2,'c','k',1)
         patchNormPlot(fliplr(Fh2),Vh2);
        
        axisGeom;axis off;
        camlight headlight;        
        drawnow;
       
        [Ff1,Vf1]=joinElementSets({F1,Fh1,fliplr(Fh2)},{V1,Vh1,Vh2});
        [Ff1,Vf1]=patchCleanUnused(Ff1,Vf1);
        [Ff1,Vf1]=mergeVertices(Ff1,Vf1);

        cFigure; hold on;

         
         gpatch(Ff1,Vf1,'w','k',1)
         patchNormPlot(Ff1,Vf1);
         

        axisGeom;axis off;
        camlight headlight;        
        drawnow;
       
        Eb=patchBoundary(Ff1,Vf1);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve1=edgeListToCurve(E1);
        indCurve1=indCurve1(1:end-1);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve2=edgeListToCurve(E1);
        indCurve2=indCurve2(1:end-1);
        
         cFigure; hold on;

         plotV(Vf1(indCurve1,:),'k.-','LineWidth',4);
         plotV(Vf1(indCurve2,:),'k.-','LineWidth',4);

         gpatch(Ff1,Vf1,'w','k',0.25)
         patchNormPlot(Ff1,Vf1);
        
        axisGeom;axis off;
        camlight headlight;        
        drawnow;
 
        logicFacesOut_implant=all(logicOut_implant(F_cylinder2),2);
        logicFacesOut_implant=triSurfLogicSharpFix(F_cylinder2,logicFacesOut_implant,3);
        F1=F_cylinder2(logicFacesOut_implant,:);
        V1=V_cylinder2;
        [F1,V1]=patchCleanUnused(F1,V1);
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [V1]=patchSmooth(F1,V1,[],cParSmooth);
        Eb=patchBoundary(F1,V1);
        
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        [~,indKeep]=max(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve3=edgeListToCurve(E1);
        indCurve3=indCurve3(1:end-1);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E1=Eb(logicKeep,:);
        indCurve4=edgeListToCurve(E1);
        indCurve4=indCurve4(1:end-1);
        
        cFigure; hold on;
        gpatch(Ff1,Vf1,'w','k',0.25);
        
        plotV(Vf1(indCurve1,:),'k.-','LineWidth',3);
        plotV(Vf1(indCurve2,:),'k.-','LineWidth',3);
        
        gpatch(F1,V1,'r','k',0.25)
        
        plotV(V1(indCurve3,:),'r.-','LineWidth',3);
        plotV(V1(indCurve4,:),'r.-','LineWidth',3);
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;
        
        D1=min(minDist(Vf1(indCurve1,:),V1(indCurve3,:)));
        D2=min(minDist(Vf1(indCurve1,:),V1(indCurve4,:)));
        
        
        if D2<D1
            [~,indMin]=minDist(Vf1(indCurve1(1),:),V1(indCurve4,:));            
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(indCurve4,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(flip(indCurve4),:)).^2,2)));
            if D2a<D1a
                indCurve4=flip(indCurve4);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh1,Vh1]=polyLoftLinear(Vf1(indCurve1,:),V1(indCurve4,:),cPar);
        else
           [~,indMin]=minDist(Vf1(indCurve1(1),:),V1(indCurve3,:));            
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(indCurve3,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve1,:)-V1(flip(indCurve3),:)).^2,2)));
            if D2a<D1a
                indCurve3=flip(indCurve3);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh1,Vh1]=polyLoftLinear(Vf1(indCurve1,:),V1(indCurve3,:),cPar);
        end
        
        D1=min(minDist(Vf1(indCurve2,:),V1(indCurve3,:)));
        D2=min(minDist(Vf1(indCurve2,:),V1(indCurve4,:)));

        if D2<D1
            [~,indMin]=minDist(Vf1(indCurve2(1),:),V1(indCurve4,:));            
            if indMin>1
                indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(indCurve4,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(flip(indCurve4),:)).^2,2)));
            if D2a<D1a
                indCurve4=flip(indCurve4);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh2,Vh2]=polyLoftLinear(Vf1(indCurve2,:),V1(indCurve4,:),cPar);
        else
           [~,indMin]=minDist(Vf1(indCurve2(1),:),V1(indCurve3,:));            
            if indMin>1
                indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
            end
            D1a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(indCurve3,:)).^2,2)));
            D2a=sum(sqrt(sum((Vf1(indCurve2,:)-V1(flip(indCurve3),:)).^2,2)));
            if D2a<D1a
                indCurve3=flip(indCurve3);
            end
            cPar.closeLoopOpt=1;
            cPar.patchType='tri_slash';
            [Fh2,Vh2]=polyLoftLinear(Vf1(indCurve2,:),V1(indCurve3,:),cPar);
        end
        
        
        cFigure; hold on;
        gpatch(Ff1,Vf1,'w','k',0.25);
        patchNormPlot(Ff1,Vf1);

        gpatch(fliplr(Fh1),Vh1,'r','k',1)
        patchNormPlot(fliplr(Fh1),Vh1);  
        
        gpatch(fliplr(Fh2),Vh2,'g','k',1)
        patchNormPlot(fliplr(Fh2),Vh2);  
        
        gpatch(F1,V1,'w','k',0.25)
        patchNormPlot(F1,V1);
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;
        
        [F_bone,V_bone]=joinElementSets({Ff1,fliplr(Fh1),fliplr(Fh2),F1},{Vf1,Vh1,Vh2,V1});
        [F_bone,V_bone]=patchCleanUnused(F_bone,V_bone);
        [F_bone,V_bone]=mergeVertices(F_bone,V_bone);
        
        cFigure; hold on;
        gpatch(F_bone,V_bone,'w','k',1);
        %patchNormPlot(F_bone,V_bone);  
        
        axisGeom;axis off;
        camlight headlight;
        drawnow;
        
        %Find intersection of cylindrical bars with muscle
        %Remove intersected points correcting the mesh. 
        % Determine optional voxel size input from mean edge size
        [D1]=patchEdgeLengths(F_muscle,V_muscle);
        [D2]=patchEdgeLengths(F_bone,V_bone);
        d=mean([D1(:);D2(:)]);
        voxelSize=d/5;
        
        % Find points outside of a simplex Cilinder1 & skin surface
        [regionLabel1]=simplexImIntersect(F_bone,V_bone,[],V_muscle,voxelSize);
        logicOut=isnan(regionLabel1);
        logicFacesOut=all(logicOut(F_muscle),2);
        logicFacesOut=triSurfLogicSharpFix(F_muscle,logicFacesOut,3);
        
        % Delet points outside of a simplex cylinder 1 and skin 
        indRim=unique(F_muscle(logicFacesOut,:));
        logicFacesSmooth=any(ismember(F_muscle,indRim),2);
        indSmooth=F_muscle(logicFacesSmooth,:);
        indRigid=F_muscle(~ismember(F_muscle,indSmooth));

        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        cParSmooth.RigidConstraints=indRigid;
        [Vm]=patchSmooth(F_muscle(logicFacesOut,:),V_muscle,[],cParSmooth);
        Fm=F_muscle(logicFacesOut,:);
        
        %%
        % Visualization
        cFigure; hold on;
        gpatch(F_bone,V_bone,'w','none',0.5);
        gpatch(Fm,Vm,'kw','k',0.5);
        %patchNormPlot(Fm,Vm);  
        axisGeom;axis off;
        camlight headlight;
        drawnow;
        
               
        %% Find edges of the muscle
        Eb=patchBoundary(Fm,Vm);
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E1=Eb(logicKeep,:);
        indRim=edgeListToCurve(E1);
        indRim=indRim(1:end-1);
        
        %Visualization
        cFigure;hold on;
        
        gpatch(F_bone,V_bone,'w','none',0.5);
        gpatch(Fm,Vm,'w','k',0.5);
        plotV(Vm(indRim,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

        
       %% Create the outer surface of the muscle by scaling the muscle
        Vma=Vm.*scaleFactor;
        Eb=patchBoundary(Fm,Vma);
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E1=Eb(logicKeep,:);
        indRima=edgeListToCurve(E1);
        indRima=indRima(1:end-1);
        
        %Visualization
        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.5);
        gpatch(Fm,Vma,'w','k',0.5);
        plotV(Vm(indRim,:),'k-','LineWidth',4);
        plotV(Vma(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;
  
        difference=mean(Vma(indRima,:),1)-mean(Vm(indRim,:),1);
        Vma=Vma-difference;
      
        %Visualization
        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.25);
        %patchNormPlot(Fm,Vm);  
        gpatch(Fm,Vma,'w','k',0.25);
        %patchNormPlot(Fm,Vma); 
        plotV(Vm(indRim,:),'k-','LineWidth',4);
        plotV(Vma(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

        [Fmb,Vmb]=regionTriMesh3D({Vm(indRim,:),Vma(indRima,:)},pointSpacing,0,'natural');
        Fmb=fliplr(Fmb);
        %Visualization
        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.25);
        %patchNormPlot(Fm,Vm);  
        gpatch(Fm,Vma,'w','k',0.25);
        %patchNormPlot(Fm,Vma); 
        gpatch(Fmb,Vmb,'w','k',0.25);
        %patchNormPlot(Fmb,Vmb); 
        
        plotV(Vm(indRim,:),'k-','LineWidth',4);
        plotV(Vma(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

     
        %% Select the most lower edge of the amputated muscle to build the
        % support structure for the mold.     
        [C,I]=min(Vma(:,3));
        snapTolerance=mean(patchEdgeLengths(Fm,Vma))/100;
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P=Vma(I,:);%Point on plane
        P(:,3)=P(:,3)+CutHeight;
        C=[];
        %Slicing surface 
        [Fmc,Vmc,~,logicSide,Eb]=triSurfSlice(Fm,Vma,C,P,n,snapTolerance);
        %Compose isolated cut geometry and boundary curves
        [Fmc,Vmc]=patchCleanUnused(Fmc(~logicSide,:),Vmc);
        
        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.25);
        patchNormPlot(Fm,Vm);  
        gpatch(Fmc,Vmc,'w','k',0.25);
        patchNormPlot(Fmc,Vmc); 
        gpatch(Fmb,Vmb,'w','k',0.25);
        patchNormPlot(Fmb,Vmb); 
        
        plotV(Vm(indRim,:),'k-','LineWidth',4);
        plotV(Vma(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;
        
        %Join and merge surfaces
        [FM,VM]=joinElementSets({Fm,Fmb,fliplr(Fmc)},{Vm,Vmb,Vmc});
        [FM,VM]=patchCleanUnused(FM,VM);
        [FM,VM]=mergeVertices(FM,VM);
        
        cFigure;hold on;
        
        gpatch(FM,VM,'w','none',1);
        patchNormPlot(FM,VM);  
        axisGeom; axis off;
        camlight headlight;
        drawnow;

        %Find edges of the outer muscle surface
        Eb=patchBoundary(FM,VM);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct); 
        
        logicKeep=G==1;
        Eb_keep=Eb(logicKeep,:);
        indRim=edgeListToCurve(Eb_keep);
        indRim=indRim(1:end-1);

        %Visualization
        cFigure;hold on;
        
        gpatch(FM,VM,'w','none',1);
        plotV(VM(indRim,:),'k-','LineWidth',4);
        axisGeom;axis off;
        camlight headlight;   
        drawnow;


        %% Create the support structure for the mold
        % The bottom curve of the outer muscle surface is scaled
        V1=VM(indRim,:).*scaleFactor1;
        difference=mean(V1,1)-mean(VM(indRim,:),1);
        V1=V1-difference;
        V1(:,3)=V1(:,3)-StandHeight;        
            
        cFigure;hold on;
        
        gpatch(FM,VM,'w','none',1);
        plotV(VM(indRim,:),'k-','LineWidth',4);
        plotV(V1,'k-','LineWidth',4);
                
        axisGeom; axis off;
        camlight headlight;   
        drawnow;

        %Mesh the space between two curves by creating the loft between
        %them
        cPar.closeLoopOpt=1;
        cPar.patchType='tri';
        [F2,V2,~,~]=polyLoftLinear(VM(indRim,:),V1,cPar);

        pointSpacing=mean(patchEdgeLengths(F2,V2));
        [F3,V3]=regionTriMesh3D({V1},pointSpacing,0,'linear');
        
        cFigure;hold on;
        
        gpatch(FM,VM,'w','none',0.5);
        %patchNormPlot(FM,VM); 
        gpatch(F3,V3,'w','k',0.5);
        %patchNormPlot(F3,V3); 
        gpatch(F2,V2,'w','k',0.5);
        %patchNormPlot(F2,V2); 
    
        plotV(VM(indRim,:),'k-','LineWidth',4);
        plotV(V1,'k-','LineWidth',4);
                
        axisGeom; axis off;
        camlight headlight;   
        drawnow;
  
        % Join and merge surfaces
        [FM,VM]=joinElementSets({FM,fliplr(F2),fliplr(F3)},{VM,V2,V3});
        [FM,VM]=patchCleanUnused(FM,VM);
        [FM,VM]=mergeVertices(FM,VM);
        
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [VM]=patchSmooth(FM,VM,[],cParSmooth);
        
        %Visualization
        cFigure;hold on;
        
        %gpatch(F_bone,V_bone,'w','none',0.25);
        gpatch(FM,VM,'w','none',0.5);        
        %patchNormPlot(FM,VM);
        
        axisGeom;axis off;
        camlight headlight;   
        drawnow;       

        %% Cut the mould into two parts for the nylon 3d printing
        snapTolerance=mean(patchEdgeLengths(FM,VM))/100; %Tolerance for surface slicing
        n=vecnormalize([1 0 0]); %Normal direction to plane
        P_cut=mean(VM,1); %Point on plane; %Point on plane
        %Slicing surface
        [Fc,Vc,~,logicSide,Eb]=triSurfSlice(FM,VM,[],P_cut,n);
        indSliceCurve=edgeListToCurve(Eb);
        indSliceCurve=indSliceCurve(1:end-1);
        Vc_slice=Vc(indSliceCurve,:);
        
        
        FM1=Fc(logicSide==1,:);
        [FM1,VM1]=patchCleanUnused(FM1,Vc);
        [FMM,VMM]=regionTriMesh3D({Vc_slice},[],0,'linear');
        FM2=Fc(logicSide==0,:);
        [FM2,VM2]=patchCleanUnused(FM2,Vc);
        
        cFigure; subplot(1,2,1);
        hold on;
        gpatch(FM1,VM1,'w','none',0.25);
        gpatch(FMM,VMM,'bw','k',1);
        patchNormPlot(fliplr(FMM),VMM);
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        
        subplot(1,2,2); hold on;
        gpatch(FM2,VM2,'w','none',0.25);
        gpatch(FMM,VMM,'bw','k',1);
        patchNormPlot(FMM,VMM);
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        drawnow;

        [FM1,VM1]=joinElementSets({FM1,fliplr(FMM)},{VM1,VMM});
        [FM1,VM1]=patchCleanUnused(FM1,VM1);
        [FM1,VM1]=mergeVertices(FM1,VM1);
        
        [FM2,VM2]=joinElementSets({FM2,FMM},{VM2,VMM});
        [FM2,VM2]=patchCleanUnused(FM2,VM2);
        [FM2,VM2]=mergeVertices(FM2,VM2);

        cFigure; hold on;
        gpatch(FM1,VM1,'w','k',1);
        patchNormPlot(FM1,VM1);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;
        
        cFigure; hold on;
        gpatch(FM2,VM2,'w','k',1);
        patchNormPlot(FM2,VM2);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;

        
        %% Save model
        if saveOn==1
            switch select_amputation_case
                case 'tf'
                    FT_mold{1}=F_bone;
                    VT_mold{1}=V_bone;
                    CT_mold{1}=1*ones(size(F_bone,1),1);
                    
                    FT_mold{2}=FM1;
                    VT_mold{2}=FM1;
                    CT_mold{2}=2*ones(size(FM1,1),1);
                    
                    FT_mold{3}=FM2;
                    VT_mold{3}=VM2;
                    CT_mold{3}=3*ones(size(FM2,1),1);
                    
                    saveName_mat=fullfile(saveFolder_mat,[saveNameGeom_tf,'.mat']);
                    save(saveName_mat,'FT_mold','VT_mold','CT_mold');
            end
        end
        
       
end
