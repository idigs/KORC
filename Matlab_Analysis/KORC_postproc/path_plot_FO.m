%% Load KORC output file

load('ST.mat')


%% Analyzing output

X=squeeze(ST.data.sp1.X(1,1,:));
Y=squeeze(ST.data.sp1.X(1,2,:));
Z=squeeze(ST.data.sp1.X(1,3,:));

if strcmp(ST.params.simulation.field_model(1),'ANALYTICAL')
    FLAGRZ=[1;1];    
else
    FLAG2D=ST.params.fields.Flag2D;
    RF=ST.params.fields.R;
    ZF=ST.params.fields.Z;

    FLAGRZ=contourc(RF,ZF,FLAG2D',1);
    FLAGRZ(:,1)=[];    
end

[PHI,R,Z]=cart2pol(X,Y,Z);

time=ST.time;

%% Plotting

orbitplots=1;
poincareplot=1;

if orbitplots==1
    
    figure;
    plot3(X,Y,Z)
    hold on; 
    plot3(X(1),Y(1),Z(1),'ko','MarkerFaceColor',[0,0,0]);
    
    plot3(FLAGRZ(1,2:end),zeros(size(FLAGRZ(1,2:end))),FLAGRZ(2,2:end))
    plot3(zeros(size(FLAGRZ(1,2:end))),FLAGRZ(1,2:end),FLAGRZ(2,2:end))
    plot3(-FLAGRZ(1,2:end),zeros(size(FLAGRZ(1,2:end))),FLAGRZ(2,2:end))
    plot3(zeros(size(FLAGRZ(1,2:end))),-FLAGRZ(1,2:end),FLAGRZ(2,2:end))
    
    hold off;
    xlabel('X')
    ylabel('y')
    zlabel('Z')
    daspect([1,1,1])
    
end

if poincareplot==1
    
    figure;
    scatter(R,Z,'.')
    hold on; 
    scatter(R(1),Z(1),'ko','MarkerFaceColor',[0,0,0]);  
    plot(FLAGRZ(1,2:end),FLAGRZ(2,2:end))
    hold off;
    xlabel('R')
    ylabel('Z')
    daspect([1,1,1])
    
end
