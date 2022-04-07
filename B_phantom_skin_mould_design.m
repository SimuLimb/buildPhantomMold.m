clear; close all; clc;

%% Build the phantom mould for muscle surface

%% Control parameters 
% Path names
projectFolder = fileparts(fileparts(mfilename('fullpath')));
projectFolder_sim_Amputation=fileparts('/Users/s2986149/Desktop/');
loadFolder=fullfile(projectFolder_sim_Amputation,'simulateAmputation.m','data','BodyParts3D','post'); 
saveFolder_stl=fullfile(projectFolder,'data','mould_stl'); 
saveFolder_mat=fullfile(projectFolder,'data','mould_mat'); 
loadGeom_tf='BodyParts3D_right_leg_transfemoral_amp_muscle_mould';
saveNameGeom_tf='BodyParts3D_right_leg_transfemoral_amp_skin_mould';

saveOn=1;
green = 1/255*[0, 100, 0];

%Geometric parameters
pointSpacing=5;
scaleFactor=[1.1 1.1 1.01];
scaleFactor1=[1.5 1.5 1];
CutHeight=50;
StandHeight=70;
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

fileName_mat=fullfile(saveFolder_mat,[loadGeom_tf,'.mat']);
model=load(fileName_mat);
F_mould=model.FT_mold; %Faces
V_mould=model.VT_mold; %Vertices
C_mould=model.CT_mold; %Color label

switch select_amputation_case
    case 'tf'     
        %Skin
        F_skin=F{2};
        V_skin=V{2};

    
        F_bone=F_mould{1};
        V_bone=V_mould{1};
        
        % Visualize surface and landmarks
        cFigure; hold on;
        gpatch(F_skin,V_skin,'w','none',0.5);
        patchNormPlot(F_skin,V_skin);
        gpatch(F_bone,V_bone,'w','none',1);
        patchNormPlot(F_bone,V_bone);
        

        axisGeom; axis off; camlight headlight;
        gdrawnow;

        % Find the intersection between cylinders and femur
        % Delete intersected points from the femur mesh
        
        %Find intersection of bone with skin
        %Remove intersected points correcting the mesh. 
        % Determine optional voxel size input from mean edge size
        [D1]=patchEdgeLengths(F_skin,V_skin);
        [D2]=patchEdgeLengths(F_bone,V_bone);
        d=mean([D1(:);D2(:)]);
        voxelSize=d/5;
        
        % Find points outside of a simplex Cilinder1 & skin surface
        [regionLabel1]=simplexImIntersect(F_bone,V_bone,[],V_skin,voxelSize);
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
        [Vs]=patchSmooth(F_skin(logicFacesOut,:),V_skin,[],cParSmooth);
        Fs=F_skin(logicFacesOut,:);
        
        %%
        % Visualization
        cFigure; hold on;
        gpatch(F_bone,V_bone,'w','none',0.5);
        gpatch(Fs,Vs,'kw','k',0.5);
        patchNormPlot(Fs,Vs);  
        axisGeom;axis off;
        camlight headlight;
        drawnow;
               
        %% Find edges of the muscle
        Eb=patchBoundary(Fs,Vs);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E1=Eb(logicKeep,:);
        indRim=edgeListToCurve(E1);
        indRim=indRim(1:end-1);
        
        %Visualization
        cFigure;hold on;
        
        gpatch(F_bone,V_bone,'w','none',0.5);
        gpatch(Fs,Vs,'w','k',0.5);
        plotV(Vs(indRim,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

        
       %% Create the outer surface of the muscle by scaling the muscle
        Vsa=Vs.*scaleFactor;
        Eb=patchBoundary(Fs,Vsa);
        [G,~,groupSize]=tesgroup(Eb,groupStruct);
        logicKeep=G==1;
        E1=Eb(logicKeep,:);
        indRima=edgeListToCurve(E1);
        indRima=indRima(1:end-1);
        
        %Visualization
        cFigure;hold on;
        
        gpatch(Fs,Vs,'w','none',0.5);
        gpatch(Fs,Vsa,'w','k',0.5);
        plotV(Vs(indRim,:),'k-','LineWidth',4);
        plotV(Vsa(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;
  
        difference=mean(Vsa(indRima,:),1)-mean(Vs(indRim,:),1);
        Vsa=Vsa-difference;
      
        %Visualization
        cFigure;hold on;
        
        gpatch(Fs,Vs,'w','none',0.25);
        %patchNormPlot(Fm,Vm);  
        gpatch(Fs,Vsa,'w','k',0.25);
        %patchNormPlot(Fm,Vma); 
        plotV(Vs(indRim,:),'k-','LineWidth',4);
        plotV(Vsa(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

        [Fsb,Vsb]=regionTriMesh3D({Vs(indRim,:),Vsa(indRima,:)},pointSpacing,0,'natural');
        Fsb=fliplr(Fsb);
        %Visualization
        cFigure;hold on;
        
        gpatch(Fs,Vs,'w','none',0.25);
        %patchNormPlot(Fm,Vm);  
        gpatch(Fs,Vsa,'w','k',0.25);
        %patchNormPlot(Fm,Vma); 
        gpatch(Fsb,Vsb,'w','k',0.25);
        %patchNormPlot(Fmb,Vmb); 
        
        plotV(Vs(indRim,:),'k-','LineWidth',4);
        plotV(Vsa(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

     
        %% Select the most lower edge of the amputated muscle to build the
        % support structure for the mold.     
        [C,I]=min(Vsa(:,3));
        snapTolerance=mean(patchEdgeLengths(Fs,Vsa))/100;
        n=vecnormalize([0 0 1]); %Normal direction to plane
        P=Vsa(I,:);%Point on plane
        P(:,3)=P(:,3)+CutHeight;
        C=[];
        %Slicing surface 
        [Fsc,Vsc,~,logicSide,Eb]=triSurfSlice(Fs,Vsa,C,P,n,snapTolerance);
        %Compose isolated cut geometry and boundary curves
        [Fsc,Vsc]=patchCleanUnused(Fsc(~logicSide,:),Vsc);
        
        cFigure;hold on;
        
        gpatch(Fs,Vs,'w','none',0.25);
        patchNormPlot(fliplr(Fs),Vs);  
        gpatch(Fsc,Vsc,'w','k',0.25);
        patchNormPlot(Fsc,Vsc); 
        gpatch(Fsb,Vsb,'w','k',0.25);
        patchNormPlot(Fsb,Vsb); 
        
        plotV(Vs(indRim,:),'k-','LineWidth',4);
        plotV(Vsa(indRima,:),'k-','LineWidth',4);

        axisGeom; axis off;
        camlight headlight;
        drawnow;

        %Join and merge surfaces
        [FS,VS]=joinElementSets({fliplr(Fs),Fsb,Fsc},{Vs,Vsb,Vsc});
        [FS,VS]=patchCleanUnused(FS,VS);
        [FS,VS]=mergeVertices(FS,VS);
        
        cFigure;hold on;
        
        gpatch(FS,VS,'w','none',1);
        patchNormPlot(FS,VS);  
        axisGeom; axis off;
        camlight headlight;
        drawnow;

        %Find edges of the outer muscle surface
        Eb=patchBoundary(FS,VS);
        groupStruct.outputType='label';
        [G,~,groupSize]=tesgroup(Eb,groupStruct); 
        
        logicKeep=G==1;
        Eb_keep=Eb(logicKeep,:);
        indRim=edgeListToCurve(Eb_keep);
        indRim=indRim(1:end-1);

        %Visualization
        cFigure;hold on;
        
        gpatch(FS,VS,'w','none',1);
        plotV(VS(indRim,:),'k-','LineWidth',4);
        axisGeom;axis off;
        camlight headlight;   
        drawnow;


        %% Create the support structure for the mold
        % The bottom curve of the outer muscle surface is scaled
        V1=VS(indRim,:).*scaleFactor1;
        difference=mean(V1,1)-mean(VS(indRim,:),1);
        V1=V1-difference;
        V1(:,3)=V1(:,3)-StandHeight;        
            
        cFigure;hold on;
        
        gpatch(FS,VS,'w','none',1);
        plotV(VS(indRim,:),'k-','LineWidth',4);
        plotV(V1,'k-','LineWidth',4);
                
        axisGeom; axis off;
        camlight headlight;   
        drawnow;

        %Mesh the space between two curves by creating the loft between
        %them
        cPar.closeLoopOpt=1;
        cPar.patchType='tri';
        [F2,V2,~,~]=polyLoftLinear(VS(indRim,:),V1,cPar);

        pointSpacing=mean(patchEdgeLengths(F2,V2));
        [F3,V3]=regionTriMesh3D({V1},pointSpacing,0,'linear');
        
        cFigure;hold on;
        
        gpatch(FS,VS,'w','none',0.5);
        patchNormPlot(FS,VS); 
        gpatch(F3,V3,'w','k',0.5);
        patchNormPlot(F3,V3); 
        gpatch(F2,V2,'w','k',0.5);
        patchNormPlot(F2,V2); 
    
        plotV(VS(indRim,:),'k-','LineWidth',4);
        plotV(V1,'k-','LineWidth',4);
                
        axisGeom; axis off;
        camlight headlight;   
        drawnow;

        % Join and merge surfaces
        [FS,VS]=joinElementSets({FS,fliplr(F2),fliplr(F3)},{VS,V2,V3});
        [FS,VS]=patchCleanUnused(FS,VS);
        [FS,VS]=mergeVertices(FS,VS);
        
        clear cParSmooth;
        cParSmooth.n=25;
        cParSmooth.Method='HC';
        [VS]=patchSmooth(FS,VS,[],cParSmooth);
        
        %Visualization
        cFigure;hold on;
        
        gpatch(F_bone,V_bone,'w','none',0.25);
        gpatch(FS,VS,'w','none',0.5);        
        %patchNormPlot(FS,VS);
        
        axisGeom;axis off;
        camlight headlight;   
        drawnow;       

        %% Cut the mould into two parts for the nylon 3d printing
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
        patchNormPlot(fliplr(FSS),VSS);
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        
        subplot(1,2,2); hold on;
        gpatch(FS2,VS2,'w','none',0.25);
        gpatch(FSS,VSS,'bw','k',1);
        patchNormPlot(FSS,VSS);
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        set(gca,'FontSize',25);
        drawnow;

        [FS1,VS1]=joinElementSets({FS1,fliplr(FSS)},{VS1,VSS});
        [FS1,VS1]=patchCleanUnused(FS1,VS1);
        [FS1,VS1]=mergeVertices(FS1,VS1);
        
        [FS2,VS2]=joinElementSets({FS2,FSS},{VS2,VSS});
        [FS2,VS2]=patchCleanUnused(FS2,VS2);
        [FS2,VS2]=mergeVertices(FS2,VS2);

        cFigure; hold on;
        gpatch(FS1,VS1,'w','k',1);
        patchNormPlot(FS1,VS1);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;
        
        cFigure; hold on;
        gpatch(FS2,VS2,'w','k',1);
        patchNormPlot(FS2,VS2);
        axisGeom; axis off; axis manual; camlight headligth;
        set(gca,'FontSize',25);
        drawnow;

        
        %% Save model
        if saveOn==1
            switch select_amputation_case
                case 'tf'
                    FT_mold{1}=FS1;
                    VT_mold{1}=FS1;
                    CT_mold{1}=1*ones(size(FS1,1),1);
                    
                    FT_mold{2}=FS2;
                    VT_mold{2}=VS2;
                    CT_mold{2}=2*ones(size(FS2,1),1);
                    
                    saveName_mat=fullfile(saveFolder_mat,[saveNameGeom_tf,'.mat']);
                    save(saveName_mat,'FT_mold','VT_mold','CT_mold');
            end
        end
        
       
end
