clear; close all; clc;

%% Control parameters 
% Path names
projectFolder = fileparts(fileparts(mfilename('fullpath')));
projectFolder_sim_Amputation=fileparts('/Users/s2986149/Desktop/');
loadFolder=fullfile(projectFolder_sim_Amputation,'simulateAmputation.m','data','BodyParts3D','post'); 
saveFolder_stl=fullfile(projectFolder,'data','mould_stl'); 
saveFolder_mat=fullfile(projectFolder,'data','mould_mat'); 
saveNameGeom_tf='BodyParts3D_right_leg_transfemoral_amp_mold';

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
        %muscle and skin+fat layers laying on them

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
        
        % Refine surface using subTriDual
        D=sqrt(sum((V_skin-[-174 -78 757]).^2,2));
        logicVertices=D<20; %Vertex logic
        logicFaces=any(logicVertices(F_skin),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(F_skin,logicFaces,3);
        [Ft,Vt,C_type,indIni]=subTriDual(F_skin,V_skin,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);

        D=sqrt(sum((Vt-[-113 -161 768]).^2,2));
        logicVertices=D<20; %Vertex logic
        logicFaces=any(logicVertices(Ft),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(Ft,logicFaces,3);
        
        [Ft,Vt,C_type,indIni]=subTriDual(Ft,Vt,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);
        
        D=sqrt(sum((Vt-[2.2 -86 756]).^2,2));
        logicVertices=D<20; %Vertex logic
        logicFaces=any(logicVertices(Ft),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(Ft,logicFaces,3);
        
        [Ft,Vt,C_type,indIni]=subTriDual(Ft,Vt,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);
        
        D=sqrt(sum((Vt-[-109 19 770]).^2,2));
        logicVertices=D<20; %Vertex logic
        logicFaces=any(logicVertices(Ft),2); %Convert to face logic
        logicFaces=triSurfLogicSharpFix(Ft,logicFaces,3);
        
        [Ft,Vt,C_type,indIni]=subTriDual(Ft,Vt,logicFaces);
        
        %Smoothen newly introduced nodes
        cPar.Method='HC'; %Smoothing method
        cPar.n=50; %Number of iterations
        cPar.RigidConstraints=indIni; %Constrained points
        [Vt]=tesSmooth(Ft,Vt,[],cPar);
        
        F_skin=Ft;
        V_skin=Vt;
        
        %Visualization
        cFigure; hold on;
        gpatch(F_cylinder,V_cylinder1,'gw','none',1);
        gpatch(fliplr(F_cylinder),V_cylinder2,'gw','none',1);
        gpatch(F_skin,V_skin,'w','k',1);
        axisGeom;
        camlight headlight;
        drawnow;

        % Determine optional voxel size input from mean edge size
        [D1]=patchEdgeLengths(F_skin,V_skin);
        [D2]=patchEdgeLengths(F_cylinder,V_cylinder1);
        d=mean([D1(:);D2(:)]);
        voxelSize=d/5;
        
        % Find points outside of a simplex Cilinder1 & skin surface
        [regionLabel1]=simplexImIntersect(F_cylinder,V_cylinder1,[],V_skin,voxelSize);
        logicOut=isnan(regionLabel1);
        logicFacesOut=all(logicOut(F_skin),2);
        logicFacesOut=triSurfLogicSharpFix(F_skin,logicFacesOut,3);
        
        % Delet points outside of a simplex cylinder 1 and skin 
        indRim=unique(F_skin(logicFacesOut,:));
        logicFacesSmooth=any(ismember(F_skin,indRim),2);
        indSmooth=F_skin(logicFacesSmooth,:);
        indRigid=F_skin(~ismember(F_skin,indSmooth));

        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        cParSmooth.RigidConstraints=indRigid;
        [Vm1]=patchSmooth(F_skin(logicFacesOut,:),V_skin,[],cParSmooth);
        Fm1=F_skin(logicFacesOut,:);

        % Find points outside of a simplex Cilinder2 & skin surface1
        [regionLabel2]=simplexImIntersect(F_cylinder,V_cylinder2,[],Vm1,voxelSize);
        logicOut=isnan(regionLabel2);
        logicFacesOut=all(logicOut(Fm1),2);
        logicFacesOut=triSurfLogicSharpFix(Fm1,logicFacesOut,3);
        
        % Delet points outside of a simplex cylinder 2 and skin 
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

        % Merge skin vertices
        [Fm,Vm]=patchCleanUnused(Fm2,Vm2);
        [Fm,Vm]=mergeVertices(Fm,Vm);
               
        %Find edges of the skin
        Eb=patchBoundary(Fm,Vm);
        groupStruct.outputType='label';
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
        
        %% Create the outer surface of the skin by scaling the skin
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
      
        cFigure;hold on;

        gpatch(F_cylinder,V_cylinder1,'g','none',0.5);
        gpatch(F_cylinder,V_cylinder2,'g','none',0.5);

        gpatch(fliplr(Fm),Vm,'kw','k',0.5);
        patchNormPlot(fliplr(Fm),Vm);
        
        gpatch(fliplr(Fm),Vm2,'bw','k',0.5);
        patchNormPlot(Fm,Vm2);
        
        plotV(Vm(indRim1,:),'g-','LineWidth',4);
        plotV(Vm2(indRim1a,:),'g-','LineWidth',4);
        
        plotV(Vm(indRim2,:),'r-','LineWidth',4);
        plotV(Vm2(indRim2a,:),'r-','LineWidth',4);
        
        plotV(Vm(indRim3,:),'m-','LineWidth',4);        
        plotV(Vm2(indRim3a,:),'m-','LineWidth',4);
        axisGeom;
        camlight headlight;
        drawnow;

        % Select the most lower edge of the amputated skin to force the
        % support structure for the mold.     
        [C,I]=min(Vm2(:,3));
        snapTolerance=mean(patchEdgeLengths(Fm,Vm2))/100;
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P=Vm2(I,:);%Point on plane
        P(:,3)=P(:,3)+60;
        C=[];
        %Slicing surface 
        [Fm3,Vm3,~,logicSide,Eb]=triSurfSlice(Fm,Vm2,C,P,n,snapTolerance);
        %Compose isolated cut geometry and boundary curves
        [Fm3,Vm3]=patchCleanUnused(Fm3(~logicSide,:),Vm3);
        
        %Visualization
        cFigure;hold on;

        gpatch(fliplr(Fm),Vm,'kw','none',0.5);
        patchNormPlot(fliplr(Fm),Vm);
        gpatch(Fm3,Vm3,'gw','k',1);
        patchNormPlot(Fm3,Vm3);
        
        axisGeom;
        camlight headlight;
        drawnow;

        %Find edges of the outer skin surface
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
        
        cFigure;hold on;
        
        gpatch(Fm,Vm,'kw','none',0.5);
        gpatch(Fm3,Vm3,'kw','none',0.5);

        
        plotV(Vm(indRim1,:),'g-','LineWidth',4);
        plotV(Vm3(indRim1b,:),'g-','LineWidth',4);
        
        plotV(Vm(indRim2,:),'r-','LineWidth',4);
        plotV(Vm3(indRim2b,:),'r-','LineWidth',4);
        
        plotV(Vm(indRim3,:),'m-','LineWidth',4);  
        plotV(Vm3(indRim3b,:),'m-','LineWidth',4);

        plotV(Vm3(indRim4b,:),'b-','LineWidth',4);

        axisGeom;
        camlight headlight;   
        drawnow;
        
        % Lofting the edges together
        pointSpacing=mean(patchEdgeLengths(Fm,Vm));
        %
        [~,indMin]=minDist(Vm(indRim1,:),Vm3(indRim1b,:)); 
        if indMin>1
            indRim1b=[indRim1b(indMin:end) indRim1b(1:indMin-1)];
        end

        cPar.closeLoopOpt=1;
        cPar.patchType='tri_slash';
        [Fb1,Vb1]=polyLoftLinear(Vm(indRim1,:),Vm3(indRim1b,:),cPar);
        %
        [~,indMin]=minDist(Vm(indRim2,:),Vm3(indRim2b,:)); 
        if indMin>1
            indRim2b=[indRim2b(indMin:end) indRim2b(1:indMin-1)];
        end
        %[Fb2,Vb2]=regionTriMesh3D({Vm(indRim2,:),Vm3(indRim2b,:)},pointSpacing,0,'natural');
        [Fb2,Vb2]=polyLoftLinear(Vm(indRim2,:),Vm3(indRim2b,:),cPar);
        
        %
        [~,indMin]=minDist(Vm(indRim3,:),Vm3(indRim3b,:)); 
        if indMin>1
            indRim3b=[indRim3b(indMin:end) indRim3b(1:indMin-1)];
        end
        %[Fb3,Vb3]=regionTriMesh3D({Vm(indRim3,:),Vm3(indRim4b,:)},pointSpacing,0,'natural');
        [Fb3,Vb3]=polyLoftLinear(Vm(indRim3,:),Vm3(indRim3b,:),cPar);

        %Visualization
        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.5);
        gpatch(Fm3,Vm3,'w','none',0.5);

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

        axisGeom;
        camlight headlight;   
        drawnow;

        %Create the support structure for the mold
        % The bottom curve of the outer skin surface is scaled
        V_bot1=Vm3(indRim4b,:).*scaleFactor1;
        difference=mean(V_bot1,1)-mean(Vm3(indRim4b,:),1);
        V_bot1=V_bot1-difference;
        V_bot1(:,3)=V_bot1(:,3)-80;        
            
        %Visualization
        cFigure;hold on;
        
        gpatch(Fm,Vm,'w','none',0.5);
        gpatch(Fm3,Vm3,'w','none',0.5);

        gpatch(Fb1,Vb1,'kw','k',0.5);
        gpatch(Fb2,Vb2,'kw','k',0.5);
        gpatch(Fb3,Vb3,'kw','k',0.5);

        plotV(Vm3(indRim4b,:),'b-','LineWidth',4);
        plotV(V_bot1,'b-','LineWidth',4);
                
        axisGeom;
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
                
        gpatch(Fm,Vm,'w','k',0.5);
        gpatch(Fm3,Vm3,'w','k',1);        
        gpatch(Fb1,Vb1,'w','k',0.5);
        gpatch(Fb2,Vb2,'w','k',0.5);
        gpatch(Fb3,Vb3,'w','k',0.5);
        gpatch(Fb4,Vb4,'w','k',0.5);
        gpatch(Fb5,Vb5,'w','k',0.5);
        
        axisGeom;
        camlight headlight;   
        drawnow;
            
        % Join and merge surfaces
        [FS,VS]=joinElementSets({fliplr(Fm),Fm3,Fb1,Fb2,Fb3,fliplr(Fb4),fliplr(Fb5)},{Vm,Vm3,Vb1,Vb2,Vb3,Vb4,Vb5});

        [FS,VS]=patchCleanUnused(FS,VS);
        [FS,VS]=mergeVertices(FS,VS);
        
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [VS]=patchSmooth(FS,VS,[],cParSmooth);
        
        %Visualization
        cFigure;hold on;
        gpatch(FS,VS,'w','k',1);        
        patchNormPlot(FS,VS);
        axisGeom;
        camlight headlight;   
        drawnow;    

        %
        snapTolerance=mean(patchEdgeLengths(FS,VS))/100; %Tolerance for surface slicing
        n=vecnormalize([1 0 0]); %Normal direction to plane
        P_cut=mean(VS,1); %Point on plane; %Point on plane
        %Slicing surface
        [Fc,Vc,~,logicSide,Eb]=triSurfSlice(FS,VS,[],P_cut,n);
        indSliceCurve=edgeListToCurve(Eb);
        indSliceCurve=indSliceCurve(1:end-1);
        Vc_slice=Vc(indSliceCurve,:);
        
        
        FS1=Fc(logicSide==1,:);
        [FS1,VS1]=patchCleanUnused(FS1,Vc);
        [FSS,VSS]=regionTriMesh3D({Vc_slice},[],0,'linear');
        FS2=Fc(logicSide==0,:);
        [FS2,VS2]=patchCleanUnused(FS2,Vc);
        
        cFigure; subplot(1,2,1);
        hold on;
        gpatch(FS1,VS1,'w','none',0.25);
        gpatch(FSS,VSS,'bw','k',1);
        patchNormPlot(FSS,VSS);
       
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        
        subplot(1,2,2); hold on;
        gpatch(FS2,VS2,'w','none',0.25);
        gpatch(FSS,VSS,'bw','k',1);
        patchNormPlot(fliplr(FSS),VSS);
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        drawnow;
        
        [FS1,VS1]=joinElementSets({FS1,FSS},{VS1,VSS});
        [FS1,VS1]=patchCleanUnused(FS1,VS1);
        [FS1,VS1]=mergeVertices(FS1,VS1);
        
        [FS2,VS2]=joinElementSets({FS2,fliplr(FSS)},{VS2,VSS});
        [FS2,VS2]=patchCleanUnused(FS2,VS2);
        [FS2,VS2]=mergeVertices(FS2,VS2);

        cFigure; hold on;
        gpatch(FS1,VS1,'w','k',1);
        %patchNormPlot(FS1,VS1);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;
        
        cFigure; hold on;
        gpatch(FS2,VS2,'w','k',1);
        %patchNormPlot(FS2,VS2);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;
        
        cFigure; hold on;
        gpatch(FS,VS,'w','k',1); 
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;

        
        %% Save model
        if saveOn==1
            switch select_amputation_case
                case 'tf'                    
                    FT_mold{1}=FS1;
                    VT_mold{1}=VS1;
                    CT_mold{1}=1*ones(size(FS1,1),1);
                    
                    FT_mold{2}=FS2;
                    VT_mold{2}=VS2;
                    CT_mold{2}=2*ones(size(FS2,1),1);
                                        
                    saveName_mat=fullfile(saveFolder_mat,[saveNameGeom_tf,'.mat']);
                    save(saveName_mat,'FT_mold','VT_mold','CT_mold');
            
            end
        end
        
       
end
