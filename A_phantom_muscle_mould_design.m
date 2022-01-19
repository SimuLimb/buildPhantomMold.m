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

%Geometric parameters
pointSpacing=5;
scaleFactor=[1.1 1.1 1.01];
scaleFactor1=[1.5 1.5 1];


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
        V_cylinder2(:,3)=V_cylinder2(:,3)-3*inputStruct.cylRadius;
        %%
        % Visualize surface and landmarks
        cFigure; hold on;
        gpatch(F_femur,V_femur,'w','k',0.5);
        gpatch(F_muscle,V_muscle,'w','k',0.5);
        gpatch(F_skin,V_skin,'w','k',0.5);
        plotV(Vc,'b-','LineWidth',4);
        plotV(mean(Vc,1),'b.','Markersize',30);
        
        gpatch(F_cylinder,V_cylinder1,'gw','g',1);
        patchNormPlot(F_cylinder,V_cylinder1);
        
        gpatch(F_cylinder,V_cylinder2,'gw','g',1);
        patchNormPlot(F_cylinder,V_cylinder2);
        
        axisGeom; camlight headlight;
        gdrawnow;

        %% Create holes within the femur to accomodate retractable cylinders
        % Find the intersection between cylinders and femur
        % Delete intersected points from the femur mesh
        
        %% First cylinder
        pointSpacing=mean(patchEdgeLengths(F_femur,V_femur));
        voxelSize=pointSpacing;
        [regionLabel_implant]=simplexImIntersect(F_femur,V_femur,[],V_cylinder1,voxelSize);
        logicOut_implant=~(~isnan(regionLabel_implant));
        logicFacesOut_implant=all(logicOut_implant(F_cylinder),2);
        
        [regionLabel_bone]=simplexImIntersect(F_cylinder,V_cylinder1,[],V_femur,voxelSize);
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
        Vc1=V1(indCurve1,:);
        
        %[~,indKeep]=min(groupSize);
        logicKeep=G~=indKeep;
        E1=Eb(logicKeep,:);
        indCurve2=edgeListToCurve(E1);
        indCurve2=indCurve2(1:end-1);
        Vc2=V1(indCurve2,:);

        % Visualization
        cFigure; hold on;
        gpatch(F1,V1,'w','k',1);
        
        plotV(Vc2,'g.-','LineWidth',1);        
        plotV(Vc1,'r.-','LineWidth',1);
%         
        axisGeom;
        camlight headlight;
        icolorbar;
        drawnow;

        logicFacesIn_implant=all(~logicOut_implant(F_cylinder),2);
        logicFacesIn_implant=triSurfLogicSharpFix(F_cylinder,logicFacesIn_implant,3);
        F2=F_cylinder(logicFacesIn_implant,:);
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
        Vc3=V2(indCurve3,:);
        
        [~,indKeep]=min(groupSize);
        logicKeep=G==indKeep;
        E2=Eb(logicKeep,:);
        indCurve4=edgeListToCurve(E2);
        indCurve4=indCurve4(1:end-1);
        Vc4=V2(indCurve4,:);
        
        %Visualization
        cFigure; hold on;
        gpatch(F1,V1,'w','k',1);
        
        plotV(V1(indCurve2,:),'g.-','LineWidth',1);
        plotV(V2(indCurve3,:),'m.-','LineWidth',1);
        
        plotV(V1(indCurve1,:),'r.-','LineWidth',1);
        plotV(V2(indCurve4,:),'c.-','LineWidth',1);
        gpatch(F2,V2,'r','k',1)
%         
        axisGeom;
        camlight headlight;
        icolorbar;
        drawnow;

        % Closing up holes between cut cylinders and holes withing the
        % femur
        [~,indMin]=minDist(V1(indCurve1,:),V2(indCurve3,:));
        
        if indMin>1
            indCurve3=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
        end
        
        [~,indMin]=minDist(V1(indCurve2,:),V2(indCurve4,:));
        
        if indMin>1
            indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
        end

        
        [Fh,Vh]=regionTriMesh3D({V1(indCurve1,:),V2(indCurve3,:)},pointSpacing,0,'natural');
        [Fh1,Vh1]=regionTriMesh3D({V1(indCurve2,:),V2(indCurve4,:)},pointSpacing,0,'natural');
        
        %Visualization
        cFigure; hold on;
        gpatch(F1,V1,'y','k',1);
        patchNormPlot(F1,V1);
        %plotV(V1(indCurve2,:),'c.-','LineWidth',1);
        %plotV(V2(indCurve3,:),'c.-','LineWidth',1);
        
        gpatch(F2,V2,'r','k',1)
        patchNormPlot(fliplr(F2),V2);
        gpatch(Fh,Vh,'w','k',1);
        patchNormPlot(fliplr(Fh),Vh);
        
        plotV(V1(indCurve2,:),'g.-','LineWidth',1);
        plotV(V2(indCurve4,:),'g.-','LineWidth',1);
        gpatch(Fh1,Vh1,'w','k',1);
        patchNormPlot(Fh1,Vh1);       

        axisGeom;
        camlight headlight;
        icolorbar;
        drawnow;

        %
        [Ff1,Vf1]=joinElementSets({F1,fliplr(F2),fliplr(Fh),Fh1},{V1,V2,Vh,Vh1});
        [Ff1,Vf1]=patchCleanUnused(Ff1,Vf1);
        [Ff1,Vf1]=mergeVertices(Ff1,Vf1);

        %Visualization
        cFigure; hold on;
        gpatch(Ff1,Vf1,'y','k',1);
        patchNormPlot(Ff1,Vf1);  
        
        %gpatch(F_cylinder,V_cylinder1,'w','k',1);
        %patchNormPlot(F_cylinder,V_cylinder1); 
        axisGeom;
        camlight headlight;
        drawnow;

        %% Second cylinder
        pointSpacing=mean(patchEdgeLengths(Ff1,Vf1));
        voxelSize=pointSpacing;
        [regionLabel_implant]=simplexImIntersect(Ff1,Vf1,[],V_cylinder2,voxelSize);
        logicOut_implant=~(~isnan(regionLabel_implant));
        logicFacesOut_implant=all(logicOut_implant(F_cylinder),2);
        
        [regionLabel_bone]=simplexImIntersect(F_cylinder,V_cylinder2,[],Vf1,voxelSize);
        logicOut_bone=isnan(regionLabel_bone);
        logicFacesOut_bone=all(logicOut_bone(Ff1),2);
        
        logicFacesOut_bone=triSurfLogicSharpFix(Ff1,logicFacesOut_bone,3);
        F3=Ff1(logicFacesOut_bone,:);
        V3=Vf1;
        [F3,V3]=patchCleanUnused(F3,V3);
        Eb=patchBoundary(F3,V3);
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

        logicFacesIn_implant=all(~logicOut_implant(F_cylinder),2);
        logicFacesIn_implant=triSurfLogicSharpFix(F_cylinder,logicFacesIn_implant,4);
        F4=F_cylinder(logicFacesIn_implant,:);
        V4=V_cylinder2;
        [F4,V4]=patchCleanUnused(F4,V4);
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [V4]=patchSmooth(F4,V4,[],cParSmooth);
        Eb=patchBoundary(F4,V4);
        
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E2=Eb(logicKeep,:);
        indCurve3=edgeListToCurve(E2);
        indCurve3=indCurve3(1:end-1);
        
        logicKeep=G==2;
        E2=Eb(logicKeep,:);
        indCurve4=edgeListToCurve(E2);
        indCurve4=indCurve4(1:end-1);
        
        %Visualization
        cFigure; hold on;
        gpatch(F3,V3,'y','k',1);
        patchNormPlot(F3,V3);
        
        plotV(V3(indCurve2,:),'g.-','LineWidth',1);
        plotV(V4(indCurve3,:),'m.-','LineWidth',1);
        
        plotV(V3(indCurve1,:),'r.-','LineWidth',1);
        plotV(V4(indCurve4,:),'c.-','LineWidth',1);
        
        gpatch(F4,V4,'r','k',1)
        patchNormPlot(F4,V4);
        axisGeom;
        camlight headlight;
        
        drawnow;

        [~,indMin]=minDist(V3(indCurve1,:),V3(indCurve4,:));
        
        if indMin>1
            indCurve4=[indCurve4(indMin:end) indCurve4(1:indMin-1)];
        end

        [~,indMin]=minDist(V3(indCurve2,:),V4(indCurve3,:));
        
        if indMin>1
            indCurve4=[indCurve3(indMin:end) indCurve3(1:indMin-1)];
        end

        
        [Fh2,Vh2]=regionTriMesh3D({V3(indCurve1,:),V4(indCurve4,:)},pointSpacing,0,'natural');
        [Fh3,Vh3]=regionTriMesh3D({V3(indCurve2,:),V4(indCurve3,:)},pointSpacing,0,'natural');
        
        %Visualization
        cFigure; hold on;
        gpatch(F3,V3,'y','k',1);
        
        plotV(V3(indCurve2,:),'c.-','LineWidth',1);
        plotV(V4(indCurve3,:),'c.-','LineWidth',1);
        
        gpatch(F4,V4,'r','k',1)
        patchNormPlot(F4,V4);
        
        gpatch(Fh2,Vh2,'w','k',1);
        patchNormPlot(Fh2,Vh2);
        
        plotV(V3(indCurve2,:),'g.-','LineWidth',1);
        plotV(V4(indCurve4,:),'g.-','LineWidth',1);
        
        gpatch(Fh3,Vh3,'w','k',1);
        patchNormPlot(fliplr(Fh3),Vh3);       

        axisGeom;
        camlight headlight;
        drawnow;

        [F_bone,V_bone]=joinElementSets({F3,F4,Fh2,fliplr(Fh3)},{V3,V4,Vh2,Vh3});
        [F_bone,V_bone]=patchCleanUnused(F_bone,V_bone);
        [F_bone,V_bone]=mergeVertices(F_bone,V_bone);
        
        cFigure; hold on;
        gpatch(F_bone,V_bone,'w','none',1);
        %patchNormPlot(F_bone,V_bone);       

        axisGeom; axis off;
        camlight headlight;
        drawnow;
        
        
        %% We created holes within the femur, now, we have to re-create cylinders with diameter smaller than the diameter of the holes
        % making sure that cylinders can go inside of the holes
        pointSpacingBeams=pointSpacing/2;
        inputStruct.cylRadius=4.5;
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
        V_cylinder2(:,3)=V_cylinder2(:,3)-3*5;
        
        cFigure; hold on;
        
        gpatch(F_cylinder,V_cylinder1,'w','none',1);
        %patchNormPlot(F_cylinder,V_cylinder1);
        gpatch(F_cylinder,V_cylinder2,'w','none',1);
        %patchNormPlot(fliplr(F_cylinder),V_cylinder2);
        
        gpatch(F_bone,V_bone,'w','k',1);
        %patchNormPlot(F_bone,V_bone);       

        axisGeom;axis off;
        camlight headlight;
        drawnow;
   
        %% Create the holes from the intersection of cylinders with muscle surface to accomodate the cylinders.
        %Local surface refinement using subTriDual
        D=sqrt(sum((V_muscle-[-150 -77 760]).^2,2));
        logicVertices=D<20; %Vertex logic
        logicFaces=any(logicVertices(F_muscle),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(F_muscle,logicFaces,3);
        
        [Ft,Vt,C_type,indIni]=subTriDual(F_muscle,V_muscle,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);

        D=sqrt(sum((Vt-[-115 -120 770]).^2,2));
        logicVertices=D<10; %Vertex logic
        logicFaces=any(logicVertices(Ft),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(Ft,logicFaces,3);
        [Ft,Vt,C_type,indIni]=subTriDual(Ft,Vt,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);
        
        
        D=sqrt(sum((Vt-[-2 -86 756]).^2,2));
        logicVertices=D<15; %Vertex logic
        logicFaces=any(logicVertices(Ft),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(Ft,logicFaces,3);
        [Ft,Vt,C_type,indIni]=subTriDual(Ft,Vt,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);
        
        
        D=sqrt(sum((Vt-[-109 -46 769]).^2,2));
        logicVertices=D<13; %Vertex logic
        logicFaces=any(logicVertices(Ft),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(Ft,logicFaces,3);
        [Ft,Vt,C_type,indIni]=subTriDual(Ft,Vt,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);
        
        
        F_muscle=Ft;
        V_muscle=Vt;
        
        cFigure; hold on;
        
        gpatch(F_cylinder,V_cylinder1,'gw','none',1);
        gpatch(F_cylinder,V_cylinder2,'gw','none',1);
        gpatch(F_muscle,V_muscle,'w','k',1);
        %gpatch(Ft,Vt,'w','k',0.6);
        %gpatch(F_muscle(logicFaces,:),V_muscle,'gw','k');
        axisGeom;
        camlight headlight;
        drawnow;
        
        % Find the intersection of cylinders with muscle layer.
        % Delete the intersected zones to form the semi-holes for
        % cylinders to lay in.
        
        % Determine optional voxel size input from mean edge size
        [D1]=patchEdgeLengths(F_muscle,V_muscle);
        [D2]=patchEdgeLengths(F_cylinder,V_cylinder1);
        d=mean([D1(:);D2(:)]);
        voxelSize=d/3;
        
        % Find points outside of a simplex Cilinder1 & Muscle surface
        [regionLabel1]=simplexImIntersect(F_cylinder,V_cylinder1,[],V_muscle,voxelSize);
        logicOut=isnan(regionLabel1);
        logicFacesOut=all(logicOut(F_muscle),2);
        logicFacesOut=triSurfLogicSharpFix(F_muscle,logicFacesOut,3);
        
        % Delet points outside of a simplex cylinder 1 and muscle 
        indRim=unique(F_muscle(logicFacesOut,:));
        logicFacesSmooth=any(ismember(F_muscle,indRim),2);
        indSmooth=F_muscle(logicFacesSmooth,:);
        indRigid=F_muscle(~ismember(F_muscle,indSmooth));

        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        cParSmooth.RigidConstraints=indRigid;
        [Vm1]=patchSmooth(F_muscle(logicFacesOut,:),V_muscle,[],cParSmooth);
        Fm1=F_muscle(logicFacesOut,:);

        % Find points outside of a simplex Cilinder2 & Muscle surface1
        [regionLabel2]=simplexImIntersect(F_cylinder,V_cylinder2,[],Vm1,voxelSize);
        logicOut=isnan(regionLabel2);
        logicFacesOut=all(logicOut(Fm1),2);
        logicFacesOut=triSurfLogicSharpFix(Fm1,logicFacesOut,3);
        
        % Delet points outside of a simplex cylinder 2 and muscle 
        indRim=unique(Fm1(logicFacesOut,:));
        logicFacesSmooth=any(ismember(Fm1,indRim),2);
        indSmooth=Fm1(logicFacesSmooth,:);
        indRigid=Fm1(~ismember(Fm1,indSmooth));

        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        cParSmooth.RigidConstraints=indRigid;
        [Vm2]=patchSmooth(Fm1(logicFacesOut,:),Vm1,[],cParSmooth);
        Fm2=Fm1(logicFacesOut,:);
        
        %%
        % Visualization
        cFigure;
        subplot(1,3,1); hold on;
        title('Intersecting sets')
        gpatch(F_cylinder,V_cylinder1,'g','none',0.5);
        gpatch(F_cylinder,V_cylinder2,'g','none',0.5);
        gpatch(Fm2,Vm2,'kw','k',0.5);
        axisGeom;
        camlight headlight;
        
        subplot(1,3,2); hold on;
        title('Point labels')
        gpatch(F_cylinder,V_cylinder1,'g','none',0.5);
        gpatch(F_cylinder,V_cylinder2,'g','none',0.5);
        gpatch(Fm2,Vm2,'kw','none',0.5);
        scatterV(Vm2,15,regionLabel2,'filled');
        colormap gjet; icolorbar;
        axisGeom;
        camlight headlight;
        
        subplot(1,3,3); hold on;
        title('Cropped geometry')
        gpatch(F_cylinder,V_cylinder1,'g','none',0.5);
        gpatch(F_cylinder,V_cylinder2,'g','none',0.5);
        gpatch(Fm2,Vm2,'kw','k',0.5);
        axisGeom;
        camlight headlight;
        drawnow;
        
        cFigure; hold on;
        gpatch(F_cylinder,V_cylinder1,'w','none',0.5);
        gpatch(F_cylinder,V_cylinder2,'w','none',0.5);
        gpatch(Fm2,Vm2,'w','k',1);
        axisGeom;axis off;
        camlight headlight;
        drawnow;

        % Merge muscle vertices
        [Fm,Vm]=patchCleanUnused(Fm2,Vm2);
        [Fm,Vm]=mergeVertices(Fm,Vm);
               
        %% Find edges of the muscle
        Eb=patchBoundary(Fm,Vm);
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E2=Eb(logicKeep,:);
        indRim1=edgeListToCurve(E2);
        indRim1=indRim1(1:end-1);
        
        logicKeep=G==2;
        E2=Eb(logicKeep,:);
        indRim2=edgeListToCurve(E2);
        indRim2=indRim2(1:end-1);
        
        logicKeep=G==3;
        E2=Eb(logicKeep,:);
        indRim3=edgeListToCurve(E2);
        indRim3=indRim3(1:end-1);
        
         
       %% Create the outer surface of the muscle by scaling the muscle
        Vm1=Vm.*scaleFactor;
        Eb=patchBoundary(Fm,Vm1);
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E2=Eb(logicKeep,:);
        indRim1a=edgeListToCurve(E2);
        indRim1a=indRim1a(1:end-1);
        
        logicKeep=G==2;
        E2=Eb(logicKeep,:);
        indRim2a=edgeListToCurve(E2);
        indRim2a=indRim2a(1:end-1);
        
        logicKeep=G==3;
        E2=Eb(logicKeep,:);
        indRim3a=edgeListToCurve(E2);
        indRim3a=indRim3a(1:end-1);
        
        difference=mean(Vm1(indRim1a,:),1)-mean(Vm(indRim1,:),1);
        Vm2=Vm1-difference;
      
        %Visualization
        cFigure;hold on;
        
        gpatch(F_bone,V_bone,'w','none',0.5);
        gpatch(F_cylinder,V_cylinder1,'w','none',0.5);
        gpatch(F_cylinder,V_cylinder2,'w','none',0.5);

        gpatch(fliplr(Fm),Vm,'rw','k',0.5);
        %patchNormPlot(fliplr(Fm),Vm);
        
        gpatch(fliplr(Fm),Vm2,'bw','k',0.5);
        %patchNormPlot(Fm,Vm2);
        
%         plotV(Vm(indRim1,:),'g-','LineWidth',4);
%         plotV(Vm2(indRim1a,:),'g-','LineWidth',4);
%         
%         plotV(Vm(indRim2,:),'r-','LineWidth',4);
%         plotV(Vm2(indRim2a,:),'r-','LineWidth',4);
%         
%         plotV(Vm(indRim3,:),'m-','LineWidth',4);        
%         plotV(Vm2(indRim3a,:),'m-','LineWidth',4);
        axisGeom; axis off;
        camlight headlight;
        drawnow;
        
        %% Select the most lower edge of the amputated muscle to build the
        % support structure for the mold.     
        [C,I]=min(Vm2(:,3));
        snapTolerance=mean(patchEdgeLengths(Fm,Vm2))/100;
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P=Vm2(I,:);%Point on plane
        P(:,3)=P(:,3)+30;
        C=[];
        %Slicing surface 
        [Fm3,Vm3,~,logicSide,Eb]=triSurfSlice(Fm,Vm2,C,P,n,snapTolerance);
        %Compose isolated cut geometry and boundary curves
        [Fm3,Vm3]=patchCleanUnused(Fm3(~logicSide,:),Vm3);
        
        %Visualization
        cFigure;hold on;
        
        %gpatch(F_bone,V_bone,'w','none',0.5);

        gpatch(fliplr(Fm),Vm,'rw','none',0.25);
        %patchNormPlot(fliplr(Fm),Vm);
        gpatch(Fm3,Vm3,'bw','k',1);
        %patchNormPlot(Fm3,Vm3);
        
        axisGeom; axis off;
        camlight headlight;
        drawnow;

        %Find edges of the outer muscle surface
        Eb=patchBoundary(Fm3,Vm3);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct); 
        
        logicKeep=G==1;
        Eb_keep=Eb(logicKeep,:);
        indRim1b=edgeListToCurve(Eb_keep);
        indRim1b=indRim1b(1:end-1);

        logicKeep=G==2;
        Eb_keep=Eb(logicKeep,:);
        indRim2b=edgeListToCurve(Eb_keep);
        indRim2b=indRim2b(1:end-1);
        
        logicKeep=G==3;
        Eb_keep=Eb(logicKeep,:);
        indRim3b=edgeListToCurve(Eb_keep);
        indRim3b=indRim3b(1:end-1);
        
        logicKeep=G==4;
        Eb_keep=Eb(logicKeep,:);
        indRim4b=edgeListToCurve(Eb_keep);
        indRim4b=indRim4b(1:end-1);
        
        %Visualization
        cFigure;hold on;
        
        %gpatch(F_bone,V_bone,'w','none',0.5);        
        gpatch(Fm,Vm,'kw','none',0.5);
        gpatch(Fm3,Vm3,'kw','none',0.5);

        
        plotV(Vm(indRim1,:),'g-','LineWidth',4);
        plotV(Vm3(indRim1b,:),'g-','LineWidth',4);
        
        plotV(Vm(indRim2,:),'r-','LineWidth',4);
        plotV(Vm3(indRim2b,:),'r-','LineWidth',4);
        
        plotV(Vm(indRim3,:),'m-','LineWidth',4);  
        plotV(Vm3(indRim3b,:),'m-','LineWidth',4);

        plotV(Vm3(indRim4b,:),'b-','LineWidth',4);

        axisGeom;axis off;
        camlight headlight;   
        drawnow;

        %Lofting edges building up the walls of the mould
        pointSpacing=mean(patchEdgeLengths(Fm,Vm));        
        [~,indMin]=minDist(Vm(indRim1,:),Vm3(indRim1b,:)); 
        if indMin>1
            indRim1b=[indRim1b(indMin:end) indRim1b(1:indMin-1)];
        end
        %[Fb1,Vb1]=regionTriMesh3D({Vm(indRim1,:),Vm3(indRim1b,:)},pointSpacing,0,'natural');
        %%
        % Create loft
        % cPar.numSteps=17;
        cPar.closeLoopOpt=1;
        cPar.patchType='tri_slash';
        [Fb1,Vb1]=polyLoftLinear(Vm(indRim1,:),Vm3(indRim1b,:),cPar);

        [~,indMin]=minDist(Vm(indRim2,:),Vm3(indRim2b,:)); 
        if indMin>1
            indRim2b=[indRim2b(indMin:end) indRim2b(1:indMin-1)];
        end
        [Fb2,Vb2]=polyLoftLinear(Vm(indRim2,:),Vm3(indRim2b,:),cPar);

        [~,indMin]=minDist(Vm(indRim3,:),Vm3(indRim3b,:)); 
        if indMin>1
            indRim3b=[indRim3b(indMin:end) indRim3b(1:indMin-1)];
        end
        [Fb3,Vb3]=polyLoftLinear(Vm(indRim3,:),Vm3(indRim3b,:),cPar);

        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.5);
        gpatch(Fm3,Vm3,'w','none',0.5);
        %gpatch(F_bone,V_bone,'w','none',0.5);        
        gpatch(Fb1,Vb1,'kw','k',1);
        gpatch(Fb2,Vb2,'kw','k',1);
        gpatch(Fb3,Vb3,'kw','k',1);
        
        plotV(Vm(indRim1,:),'g.-','LineWidth',1);
        plotV(Vm3(indRim1b,:),'g.-','LineWidth',1);
        
        plotV(Vm(indRim2,:),'r-','LineWidth',1);
        plotV(Vm3(indRim2b,:),'r-','LineWidth',1);
        
        plotV(Vm(indRim3,:),'m-','LineWidth',1);  
        plotV(Vm3(indRim3b,:),'m-','LineWidth',1);

        plotV(Vm3(indRim4b,:),'b-','LineWidth',4);
        axisGeom; axis off;
        camlight headlight;   
        drawnow;

        %% Create the support structure for the mold
        % The bottom curve of the outer muscle surface is scaled
        V_bot1=Vm3(indRim4b,:).*scaleFactor1;
        difference=mean(V_bot1,1)-mean(Vm3(indRim4b,:),1);
        V_bot1=V_bot1-difference;
        V_bot1(:,3)=V_bot1(:,3)-50;        
            
        cFigure;hold on;
        
        %gpatch(F_bone,V_bone,'w','none',0.5);
        gpatch(Fm,Vm,'kw','none',0.5);
        gpatch(Fm3,Vm3,'kw','none',0.5);

        gpatch(Fb1,Vb1,'kw','k',0.5);
        gpatch(Fb2,Vb2,'kw','k',0.5);
        gpatch(Fb3,Vb3,'kw','k',0.5);

        plotV(Vm3(indRim4b,:),'b-','LineWidth',4);
        plotV(V_bot1,'b-','LineWidth',4);
                
        axisGeom; axis off;
        camlight headlight;   
        drawnow;

        %Mesh the space between two curves by creating the loft between
        %them
        cPar.closeLoopOpt=1;
        cPar.patchType='tri';
        [Fb4,Vb4,~,~]=polyLoftLinear(Vm3(indRim4b,:),V_bot1,cPar);

        pointSpacing=mean(patchEdgeLengths(Fb4,Vb4));
        [Fb5,Vb5]=regionTriMesh3D({V_bot1},pointSpacing,0,'linear');
        
        %Visualization
        cFigure;hold on;
        
        %gpatch(F_bone,V_bone,'w','none',0.5);
        
        gpatch(Fm,Vm,'w','none',0.5);
        gpatch(Fm3,Vm3,'w','none',0.5);        
        gpatch(Fb1,Vb1,'kw','k',0.5);
        gpatch(Fb2,Vb2,'kw','k',0.5);
        gpatch(Fb3,Vb3,'kw','k',0.5);
        gpatch(Fb4,Vb4,'kw','k',0.5);
        gpatch(Fb5,Vb5,'kw','k',0.5);
        
        plotV(Vm3(indRim4b,:),'b-','LineWidth',4);
        plotV(V_bot1,'b-','LineWidth',4);
        
        axisGeom; axis off;
        camlight headlight;   
        drawnow;

        % Join and merge surfaces
        [FM,VM]=joinElementSets({fliplr(Fm),Fm3,Fb1,Fb2,Fb3,fliplr(Fb4),fliplr(Fb5)},{Vm,Vm3,Vb1,Vb2,Vb3,Vb4,Vb5});
        [FM,VM]=patchCleanUnused(FM,VM);
        [FM,VM]=mergeVertices(FM,VM);
        
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [VM]=patchSmooth(FM,VM,[],cParSmooth);
        
        %Visualization
        cFigure;hold on;
        
        %gpatch(F_bone,V_bone,'w','none',0.5);        
        gpatch(FM,VM,'w','none',0.5);        
        %patchNormPlot(FM,VM);
        
        axisGeom;
        camlight headlight;   
        drawnow;       
 
        %Visualization
        cFigure; hold on;
        gpatch(FM,VM,'w','k',1); 
        patchNormPlot(FM,VM);
        gpatch(F_bone,V_bone,'w','k',1);
        %patchNormPlot(F_bone,V_bone);
        gpatch(F_cylinder,V_cylinder1,'w','k',1);
        patchNormPlot(F_cylinder,V_cylinder1);
        gpatch(F_cylinder,V_cylinder2,'w','k',1);
        patchNormPlot(F_cylinder,V_cylinder2);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
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
        patchNormPlot(FMM,VMM);
        
        
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        
        subplot(1,2,2); hold on;
        gpatch(FM2,VM2,'w','none',0.25);
        gpatch(FMM,VMM,'bw','k',1);
        patchNormPlot(fliplr(FMM),VMM);
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        drawnow;
        
        [FM1,VM1]=joinElementSets({FM1,FMM},{VM1,VMM});
        [FM1,VM1]=patchCleanUnused(FM1,VM1);
        [FM1,VM1]=mergeVertices(FM1,VM1);
        
        [FM2,VM2]=joinElementSets({FM2,fliplr(FMM)},{VM2,VMM});
        [FM2,VM2]=patchCleanUnused(FM2,VM2);
        [FM2,VM2]=mergeVertices(FM2,VM2);

        cFigure; hold on;
        gpatch(FM1,VM1,'w','k',1);
        %patchNormPlot(FM1,VM1);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;
        
        cFigure; hold on;
        gpatch(FM2,VM2,'w','k',1);
        %patchNormPlot(FM2,VM2);
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
                    
                    FT_mold{4}=F_cylinder;
                    VT_mold{4}=V_cylinder1;
                    CT_mold{4}=4*ones(size(F_cylinder,1),1);
                    
                    FT_mold{5}=fliplr(F_cylinder);
                    VT_mold{5}=V_cylinder2;
                    CT_mold{5}=5*ones(size(F_cylinder,1),1);

                    saveName_mat=fullfile(saveFolder_mat,[saveNameGeom_tf,'.mat']);
                    save(saveName_mat,'FT_mold','VT_mold','CT_mold');
                    
%                     stlStruct.solidNames={'Muslce_phantom'};
%                     stlStruct.solidVertices={VM};
%                     stlStruct.solidFaces={FM};
%                     stlStruct.solidNormals={[]};
%                     fileName=fullfile(saveFolder_stl,'phantom_muscle.stl');
%                     export_STL_txt(fileName,stlStruct);
            
            end
        end
        
       
end
