% Copyright : This code is written by Gerasimos Arvanitis {arvanitis@ece.upatras.gr}
%              Electrical and Computer Engineering, University of Patras. The code
%              may be used, modified and distributed for research purposes with
%              acknowledgement of the author and inclusion of this copyright information.


function handles = GUI_features()

handles.ok = 0;

% The general panel of the GUI
plot_figure.fi = figure('units','pixels',...
              'position',[50 50 1300 770],...
              'menubar','none',...
              'resize','off',...
              'numbertitle','off',...
              'name','Spectral Features Identification');

% Figure that displays the original point cloud
ax = axes('units','pixels',...
            'position',[30 200 400 400],...
            'drawmode','fast');

% Figure that displays the point cloud with the identified features         
ax1 = axes('units','pixels',...
            'position',[460 200 400 400],...
            'drawmode','fast');
        
% Figure that displays the point cloud with the identified features         
ax2 = axes('units','pixels',...
            'position',[890 200 400 400],...
            'drawmode','fast');        
        
% Editable area for inserting the preferable value of the patch size (default 15)
text = uicontrol('style','edit',...
                 'units','pix',...
                 'position',[1150 655 100 30],...
                 'fontsize',16,'string',15);

% Static text 'Size of points:'            
text1 = uicontrol('style','text',...
                 'units','pix',...
                 'position',[15 50 270 30],...
                 'fontsize',16,'string','Size of points:');   
 
% Static text 'Size of patches'              
text2 = uicontrol('style','text',...
                 'units','pix',...
                 'position',[1135 690 130 20],...
                 'fontsize',13,'string','Size of patches');              

% Text message that show if a process running or not                 
message = uicontrol('Style','text',...
         'String','   ',...
         'fontsize',19, ...
         'Position',[500 700 250 60]);

% Button for loading a new .ply file and estimating its features   
button_to_load = uicontrol('style','push',...
                 'units','pix',...
                 'position',[500 10 200 30],...
                 'backgroundcolor','w',...
                 'HorizontalAlign','left',...
                 'string','Load .ply file',...
                 'fontsize',14,'fontweight','bold',...
                 'callback',{@load_ply});

% Button for estimating the new features of a loaded .ply file based on the new selected patch size             
button_to_reload = uicontrol('style','push',...
                 'units','pix',...
                 'position',[1105 715 180 40],...
                 'backgroundcolor','w',...
                 'HorizontalAlign','left',...
                 'string','Update features',...
                 'fontsize',14,'fontweight','bold',...
                 'callback',{@pb_call_reload});  

% Checkbox for enabling or disabling the features curves displaying functionality. The default value is disabled.             
checkbox = uicontrol('Style','checkbox', ...
                     'Value',0, ...
                     'position',[1030 140 180 40], ...
                     'string','Show Feature Curves');
             
% The slider that changes the size of the points             
Slider.sl  = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 10 270 30],...
                 'min',1,'max',20,'val',4); 

% Displays the min  
Slider.ed(1) = uicontrol('style','text',...
                    'unit','pix',...
                    'position',[10 10 40 30],...
                    'fontsize',12,...
                    'string','1');   
                
% Displays the value                
Slider.ed(2) = uicontrol('style','text',...
                    'unit','pix',...
                    'position',[220 50 50 30],...
                    'fontsize',16,...
                    'string','4');  

% Displays the max                
Slider.ed(3) = uicontrol('style','text',...
                    'unit','pix',...
                    'position',[340 10 40 30],...
                    'fontsize',12,...
                    'string','20');      
% Shared Callback                
set([Slider.ed(:);Slider.sl ],'call',{@slider_change_size,Slider}); 

% List of information about the mesh
informations = uicontrol('style','text',...
                 'units','pix',...
                 'position',[950 10 320 110],...
                 'backgroundcolor','w',...
                 'string','Mesh Information',...
                 'fontsize',12,...
                 'HorizontalAlign','center');

% Button that saves the outputs (in menu area) 
button_to_save = uimenu('label','Save outputs');
set(button_to_save,'callback',{@save_results})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Function for loading a ply file and estimating features based on the default patches size %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
function [] = load_ply(varargin)

    % if there is already a mesh then clear everthing and load a new one
    if handles.ok == 1 
       message.String = 'In processing !!!'; 
       clearvars -except plot_figure.fi text Slider ax ax1 ax2 informations message checkbox 
       cla(ax)
       cla(ax1)
       cla(ax2)
       message.String = ' ';
    end
    
    handles.ok = 1;
    
    % Load only .ply files  
    [myfilename,path,indx] = uigetfile({'*.ply'});
    fullpath = strcat(path,myfilename);

    % Store the full path and the file name for using it later
    handles.fullpath = fullpath;
    handles.myfilename = myfilename;

    newvalue = get(Slider.ed(2),{'String'}); 

    [vertices,faces] = read_ply(fullpath);
    handles.vertices = vertices;
    handles.faces = double(faces);

    message.String = 'In processing !!!';

    % Display the point cloud in the figure 1
    axes(ax);
    ptCloud = pcread(fullpath);
    handles.Original(:,1) = ptCloud.Location(:,1);
    handles.Original(:,2) = ptCloud.Location(:,2);
    handles.Original(:,3) = ptCloud.Location(:,3);

    scatter3(handles.Original(:,1), handles.Original(:,2), handles.Original(:,3) , str2double(newvalue), 'b', 'filled')
    hold on
    rotate3d on
    axis tight
    title('Original');

    apoints = [ handles.vertices(:,1), handles.vertices(:,2), handles.vertices(:,3)]; 

    % Find the k nearest neighbors of each point
    anumNeighbours = str2double(get(text, 'String'));   %how many neighbors we get at every point
    apoints = double(apoints);
    akdtreeobj = KDTreeSearcher(apoints,'distance','euclidean');
    aAdjOfPoints = knnsearch(akdtreeobj,apoints,'k',(anumNeighbours+1));
    AdjOfPointsa3 = knnsearch(akdtreeobj,apoints,'k',(anumNeighbours+1)); 
    aAdjOfPoints = aAdjOfPoints(:,2:end);

    % Find the normals of each point
    [Nx,Ny,Nz]=patchnormals_double(handles.faces(:,1),handles.faces(:,2),handles.faces(:,3),double(handles.vertices(:,1)),double(handles.vertices(:,2)),double(handles.vertices(:,3)));
    handles.normalvectors=zeros(length(Nx),3);
    handles.normalvectors(:,1)=Nx;
    handles.normalvectors(:,2)=Ny;
    handles.normalvectors(:,3)=Nz;

    featuresize = size(handles.vertices,1);

    for i = 1:featuresize %for each point of the mesh

        % Create a matrix nn1 (k x 3) consisting of the normals of the k nearest neighbors of i point
        for j = 1:anumNeighbours+1
           nn1(j,:) = handles.normalvectors(AdjOfPointsa3(i,j),:);
        end

        % Estimate the matrix convn1 (3 x 3)
        convn1 = nn1'*nn1;

        % Estimate the eigenvalues of the  matrix
        [vb lb] = eig(convn1);

        calll(i,1) = lb(1,1);
        calll(i,2) = lb(2,2);
        calll(i,3) = lb(3,3);

        % Estimate the saliency of each point
        heatmap(i,1) = 1./norm(calll(i,:));
    end

    % Normalize the saliency of each point taking values [0-1]
    handles.norm_heatmap = (heatmap - min(heatmap)) / ( max(heatmap) - min(heatmap) );

    % We use kmeans for separating the normalized values in 5 clusters
    [handles.in inx] = kmedoids(handles.norm_heatmap,5); 

    % We sort the results because kmeans gives random values. However, we know that lowest values must be in class 1 and highest values in class 5
    [sort_max_inx handles.sort_max_inx_ind] = sort(inx);

    % We assume that the first two classes consists of non - features, while the three next classes consists of features with different salience 
    j =1;
    k =1;
    for i = 1:size(handles.vertices,1)
        if  handles.in(i) == handles.sort_max_inx_ind(1) || handles.in(i) == handles.sort_max_inx_ind(2)
            % The points lie in "flat" areas are represented in blue color  
            handles.blue(j,1) = handles.vertices(i,1);
            handles.blue(j,2) = handles.vertices(i,2);
            handles.blue(j,3) = handles.vertices(i,3);
            j=j+1;
        else
            % The points lie in "edges, corners, curves" areas are represented in red color  
            handles.red(k,1) = handles.vertices(i,1);
            handles.red(k,2) = handles.vertices(i,2);
            handles.red(k,3) = handles.vertices(i,3);
            red_features(k,1) = i;
            k=k+1;
        end
    end

    % Display in figure 2 the results in two different colors (i) blue (non-features), (ii) red (features)
    axes(ax1)
    scatter3(handles.blue(:,1), handles.blue(:,2), handles.blue(:,3), str2double(newvalue), 'b', 'filled')
    hold on
    scatter3(handles.red(:,1), handles.red(:,2), handles.red(:,3), str2double(newvalue), 'r', 'filled')
    hold on
    rotate3d on
    axis tight
    title('Features');
    
   %Estimate and display the feature curves only if the user has selected it. This process is time-consuming. 
   if checkbox.Value == 1
       
    %Estimate the mean curvature
    [GC MC]= curvatures(handles.vertices(:,1),handles.vertices(:,2),handles.vertices(:,3),handles.faces);
    MCMC = MC(red_features(:,1),1);
        
    %Normalize the mean curvature values
    norm_MC = (MCMC - min(MCMC)) / ( max(MCMC) - min(MCMC) );
    norm_all =  abs(norm_MC-1); 
    
    %Estimate the coefficients of the histogram of the normalized mean curvature
    [handles.coff_curv bins_curv] = hist(norm_all);
    sal = handles.norm_heatmap(red_features(:,1),1);
    [handles.coff_sal bins_sal] = hist(sal);

    % We use only these vertices which have been indicated as features
    handles.newvertices = handles.vertices(red_features(:,1),:);
    
    %Evaluate the optimal number of clusters
    clustev = evalclusters(norm_all, 'kmeans', 'CalinskiHarabasz', 'KList', 1:5);
    handles.kBest = clustev.OptimalK;
    
    %Clustering based on k-means algorithm
    [ii, iix] = kmedoids(norm_all, handles.kBest);

    %Remove small classes
    kbst = handles.kBest;
    j = 1;
       for i = 1:handles.kBest
         a(i,1) =  size(find(ii==i),1);
          if a(i,1) > floor(0.15*size(norm_all,1))
            handles.class{j} = find(ii==i);
            j = j + 1;
          else
            kbst = kbst - 1;
          end
       end
    handles.kBest = kbst;
    
    newvalue = get(Slider.ed(2),{'String'});

    %Display the results to the figure 3
    axes(ax2)
        if handles.kBest==1
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
        end
        if handles.kBest==2
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
        end
        if handles.kBest==3
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on
        end
        if handles.kBest==4
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{4},1), handles.newvertices(handles.class{4},2), handles.newvertices(handles.class{4},3), str2double(newvalue), 'y', 'filled')
            hold on
        end
        if handles.kBest==5
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{4},1), handles.newvertices(handles.class{4},2), handles.newvertices(handles.class{4},3), str2double(newvalue), 'y', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{5},1), handles.newvertices(handles.class{5},2), handles.newvertices(handles.class{5},3), str2double(newvalue), 'c', 'filled')
            hold on
        end
        rotate3d on
        axis tight
        title('Feature Curve');
    end
    
    message.String = ' ';

    set(informations,'string',{'Mesh Information'; ['Name of model: ', myfilename]; ['Vertices: ', num2str(size(handles.vertices,1))]; ['Faces: ', num2str(size(handles.faces,1))]; ['Number of Features: ', num2str(size(handles.red,1))];}) 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Function that estimates the features based on the new patches size and updates the figure %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = pb_call_reload(varargin)

    if handles.ok == 0
        
        gui_pop()

    else

        cla(ax)
        cla(ax1)
        cla(ax2)

        % The processing is starting
        message.String = 'In processing !!!';

        fullpath = handles.fullpath;
        [vertices1,faces1] = read_ply(fullpath);
        handles.vertices = vertices1;
        apoints = [ handles.vertices(:,1), handles.vertices(:,2), handles.vertices(:,3)]; 

        % We follow the same steps for the estimation of the features, as in the previous function, but we update only the figure 2
        anumNeighbours = str2double(get(text, 'String')); %how many neighbors we get at every point
        apoints = double(apoints);
        akdtreeobj = KDTreeSearcher(apoints,'distance','euclidean');
        aAdjOfPoints = knnsearch(akdtreeobj,apoints,'k',(anumNeighbours+1));
        AdjOfPointsa3 = knnsearch(akdtreeobj,apoints,'k',(anumNeighbours+1)); 
        aAdjOfPoints = aAdjOfPoints(:,2:end);

        [Nx,Ny,Nz]=patchnormals_double(handles.faces(:,1),handles.faces(:,2),handles.faces(:,3),double(handles.vertices(:,1)),double(handles.vertices(:,2)),double(handles.vertices(:,3)));
        handles.normalvectors=zeros(length(Nx),3);
        handles.normalvectors(:,1)=Nx;
        handles.normalvectors(:,2)=Ny;
        handles.normalvectors(:,3)=Nz;

        featuresize = size(handles.vertices,1);

        for i = 1:featuresize

            for j = 1:anumNeighbours+1
                nn1(j,:) = handles.normalvectors(AdjOfPointsa3(i,j),:);
            end

            convn1 = nn1'*nn1;

            [vb lb] = eig(convn1);

            calll(i,1) = lb(1,1);
            calll(i,2) = lb(2,2);
            calll(i,3) = lb(3,3);

            heatmap(i,1) = 1./norm(calll(i,:));
        end

        handles.norm_heatmap = (heatmap - min(heatmap)) / ( max(heatmap) - min(heatmap) );

        [handles.in inx] = kmedoids(handles.norm_heatmap,5); 

        [sort_max_inx handles.sort_max_inx_ind] = sort(inx);

        % There are already values in these variables so we clear them
        handles.blue = [];
        handles.red = [];

        j = 1;
        k = 1;
        for i = 1:size(handles.vertices,1)
            if  handles.in(i) == handles.sort_max_inx_ind(1) || handles.in(i) == handles.sort_max_inx_ind(2)
                handles.blue(j,1) = handles.vertices(i,1);
                handles.blue(j,2) = handles.vertices(i,2);
                handles.blue(j,3) = handles.vertices(i,3);
                j=j+1;
            else
                handles.red(k,1) = handles.vertices(i,1);
                handles.red(k,2) = handles.vertices(i,2);
                handles.red(k,3) = handles.vertices(i,3);
                red_features(k,1) = i;
                k=k+1;
            end
        end

        newvalue = get(Slider.ed(2),{'String'}); 

        axes(ax)
        scatter3(handles.Original(:,1), handles.Original(:,2), handles.Original(:,3), str2double(newvalue), 'b', 'filled')
        rotate3d on
        axis tight

        axes(ax1)
        scatter3(handles.blue(:,1), handles.blue(:,2), handles.blue(:,3), str2double(newvalue), 'b', 'filled')
        hold on 
        scatter3(handles.red(:,1), handles.red(:,2), handles.red(:,3), str2double(newvalue), 'r', 'filled')
        hold on
        rotate3d on
        axis tight
       
        
       if checkbox.Value == 1 
       [GC MC]= curvatures(handles.vertices(:,1),handles.vertices(:,2),handles.vertices(:,3),handles.faces);
        
       MCMC = MC(red_features(:,1),1);
        
       norm_MC = (MCMC - min(MCMC)) / ( max(MCMC) - min(MCMC) );
       norm_all =  abs(norm_MC-1); 
           
       [handles.coff_curv bins_curv] = hist(norm_all);
        sal = handles.norm_heatmap(red_features(:,1),1);
       [handles.coff_sal bins_sal] = hist(sal);
   
       handles.newvertices = handles.vertices(red_features(:,1),:);

       clustev = evalclusters(norm_all, 'kmeans', 'CalinskiHarabasz', 'KList', 1:5);
       handles.kBest = clustev.OptimalK;
       [ii, iix] = kmedoids(norm_all, handles.kBest);

       kbst = handles.kBest;
       j = 1;
       for i = 1:handles.kBest
         a(i,1) =  size(find(ii==i),1);
          if a(i,1) > floor(0.15*size(norm_all,1))
            handles.class{j} = find(ii==i);
            j = j + 1;
          else
            kbst = kbst - 1;
          end
       end
       handles.kBest = kbst;
       
       newvalue = get(Slider.ed(2),{'String'});

        axes(ax2)
         
        if handles.kBest==1
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
        end
        if handles.kBest==2
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
        end
        if handles.kBest==3
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on
        end
        if handles.kBest==4
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{4},1), handles.newvertices(handles.class{4},2), handles.newvertices(handles.class{4},3), str2double(newvalue), 'y', 'filled')
            hold on
        end
        if handles.kBest==5
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{4},1), handles.newvertices(handles.class{4},2), handles.newvertices(handles.class{4},3), str2double(newvalue), 'y', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{5},1), handles.newvertices(handles.class{5},2), handles.newvertices(handles.class{5},3), str2double(newvalue), 'c', 'filled')
            hold on
        end
        rotate3d on
        axis tight
        title('Feature Curve');
       end
       
        set(informations,'string',{'Mesh Information'; ['Name of model: ', handles.myfilename]; ['Vertices: ', num2str(size(handles.vertices,1))]; ['Faces: ', num2str(size(handles.faces,1))]; ['Number of Features: ', num2str(size(handles.red,1))];}) 

        message.String = " ";
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function for saving the results in output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = save_results(varargin)
    
    if  handles.ok == 0
        
        gui_pop()
        
    else  
        message.String = 'In processing !!!'; 
        
        B = regexp(handles.myfilename,'\d*','Match');
        
        % Prepare the names of the outputs
        ending = '_features.obj';
        ending1 = '_.txt';
        ending2 = '_feature_curves.obj';
        ending3 = '_coefficient.csv';
        ending4 = '_saliency_features.obj';

        
        newmesh = strcat(B{1},ending);
        newmesh1 = strcat(B{1},ending1);
        newmesh2 = strcat(B{1},ending2);
        newmesh3 = strcat(B{1},ending3);
        newmesh4 = strcat(B{1},ending4);

        % Create a .obj file showing the features in red color
        k = 1;
        fileID = fopen(newmesh, 'w');
            for i = 1:size(handles.vertices,1)
                if  handles.in(i) == handles.sort_max_inx_ind(1) || handles.in(i) == handles.sort_max_inx_ind(2)
                    fprintf(fileID,'v %f %f %f %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 0, 0, 1 );
                else
                    fprintf(fileID,'v %d %d %d %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 1, 0, 0 );  
                    red_features(k,1) = i;
                    k = k + 1;
                end
            end
            for j = 1:size(handles.faces,1)
                fprintf(fileID,'f %d %d %d\n',handles.faces(j,1), handles.faces(j,2), handles.faces(j,3)); 
            end
        fclose(fileID);
        
        % Create a .txt file presenting the indexes of the feature points
        fileID = fopen(newmesh1, 'w');
            for o = 1:size(red_features,1)
                fprintf(fileID,'%d \t',red_features(o,1));  
            end
        fclose(fileID);
        
        % Create a .obj file showing in different colors, different type of features      
        fileID = fopen(newmesh4, 'w');
            for i = 1:size(handles.vertices,1)
                if  handles.in(i) == handles.sort_max_inx_ind(1)
                    fprintf(fileID,'v %f %f %f %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 0, 0, 1 );
                elseif  handles.in(i) == handles.sort_max_inx_ind(2)
                    fprintf(fileID,'v %f %f %f %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 0, 1, 0 );
                elseif  handles.in(i) == handles.sort_max_inx_ind(3)
                    fprintf(fileID,'v %f %f %f %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 0, 1, 1 );
                elseif handles.in(i) == handles.sort_max_inx_ind(4)
                    fprintf(fileID,'v %f %f %f %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 1, 0, 1 );          
                else 
                    fprintf(fileID,'v %f %f %f %d %d %d\n',handles.vertices(i,1), handles.vertices(i,2), handles.vertices(i,3), 1, 0, 0 );         
                end
            end
            for j = 1:size(handles.faces,1)
                    fprintf(fileID,'f %d %d %d\n',handles.faces(j,1), handles.faces(j,2), handles.faces(j,3)); 
            end
        fclose(fileID);
        
        
         if checkbox.Value == 1 
         
         % Create a .csv consisting of the histogram coefficients
         coffs = [handles.coff_curv  handles.coff_sal];
         csvwrite(newmesh3,coffs);
         
         % Create a .obj file and k .txt files presenting the features curve
          if  handles.kBest==1
              
              fileID = fopen(newmesh2, 'w');
                  for i = 1:size(handles.class{1},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{1}(i,1),1), handles.newvertices(handles.class{1}(i,1),2), handles.newvertices(handles.class{1}(i,1),3), 0, 0, 1);  
                  end
              fclose(fileID); 
              
              newsubmesh1 = strcat(B{1},'_1.txt');
              fileID = fopen(newsubmesh1, 'w');
                for i = 1:size(handles.class{1},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{1}(i,1),1));  
                end
              fclose(fileID);  
              
          elseif handles.kBest==2
              
              fileID = fopen(newmesh2, 'w');
                  for i = 1:size(handles.class{1},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{1}(i,1),1), handles.newvertices(handles.class{1}(i,1),2), handles.newvertices(handles.class{1}(i,1),3), 0, 0, 1);  
                  end
                  for i = 1:size(handles.class{2},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{2}(i,1),1), handles.newvertices(handles.class{2}(i,1),2), handles.newvertices(handles.class{2}(i,1),3), 1, 0, 0);  
                  end
              fclose(fileID); 
              
              newsubmesh1 = strcat(B{1},'_1.txt');
              fileID = fopen(newsubmesh1, 'w');
                for i = 1:size(handles.class{1},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{1}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh2 = strcat(B{1},'_2.txt');
              fileID = fopen(newsubmesh2, 'w');
                for i = 1:size(handles.class{2},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{2}(i,1),1));  
                end
              fclose(fileID);
              
              
          elseif handles.kBest==3
              
              fileID = fopen(newmesh2, 'w');
                  for i = 1:size(handles.class{1},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{1}(i,1),1), handles.newvertices(handles.class{1}(i,1),2), handles.newvertices(handles.class{1}(i,1),3), 0, 0, 1);  
                  end
                  for i = 1:size(handles.class{2},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{2}(i,1),1), handles.newvertices(handles.class{2}(i,1),2), handles.newvertices(handles.class{2}(i,1),3), 1, 0, 0);  
                  end
                  for i = 1:size(handles.class{3},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{3}(i,1),1), handles.newvertices(handles.class{3}(i,1),2), handles.newvertices(handles.class{3}(i,1),3), 0, 1, 0);  
                  end
              fclose(fileID);  
              
               newsubmesh1 = strcat(B{1},'_1.txt');
              fileID = fopen(newsubmesh1, 'w');
                for i = 1:size(handles.class{1},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{1}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh2 = strcat(B{1},'_2.txt');
              fileID = fopen(newsubmesh2, 'w');
                for i = 1:size(handles.class{2},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{2}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh3 = strcat(B{1},'_3.txt');
              fileID = fopen(newsubmesh3, 'w');
                for i = 1:size(handles.class{3},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{3}(i,1),1));  
                end
              fclose(fileID);              
              
          elseif handles.kBest==4
              
              fileID = fopen(newmesh2, 'w');
                  for i = 1:size(handles.class{1},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{1}(i,1),1), handles.newvertices(handles.class{1}(i,1),2), handles.newvertices(handles.class{1}(i,1),3), 0, 0, 1);  
                  end
                  for i = 1:size(handles.class{2},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{2}(i,1),1), handles.newvertices(handles.class{2}(i,1),2), handles.newvertices(handles.class{2}(i,1),3), 1, 0, 0);  
                  end
                  for i = 1:size(handles.class{3},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{3}(i,1),1), handles.newvertices(handles.class{3}(i,1),2), handles.newvertices(handles.class{3}(i,1),3), 0, 1, 0);  
                  end
                  for i = 1:size(handles.class{4},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{4}(i,1),1), handles.newvertices(handles.class{4}(i,1),2), handles.newvertices(handles.class{4}(i,1),3), 1, 1, 0);  
                  end                  
              fclose(fileID);   
              
               newsubmesh1 = strcat(B{1},'_1.txt');
              fileID = fopen(newsubmesh1, 'w');
                for i = 1:size(handles.class{1},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{1}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh2 = strcat(B{1},'_2.txt');
              fileID = fopen(newsubmesh2, 'w');
                for i = 1:size(handles.class{2},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{2}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh3 = strcat(B{1},'_3.txt');
              fileID = fopen(newsubmesh3, 'w');
                for i = 1:size(handles.class{3},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{3}(i,1),1));  
                end
              fclose(fileID);               

              newsubmesh4 = strcat(B{1},'_4.txt');
              fileID = fopen(newsubmesh4, 'w');
                for i = 1:size(handles.class{4},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{4}(i,1),1));  
                end
              fclose(fileID); 
              
           elseif handles.kBest==5
               
              fileID = fopen(newmesh2, 'w');
                  for i = 1:size(handles.class{1},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{1}(i,1),1), handles.newvertices(handles.class{1}(i,1),2), handles.newvertices(handles.class{1}(i,1),3), 0, 0, 1);  
                  end
                  for i = 1:size(handles.class{2},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{2}(i,1),1), handles.newvertices(handles.class{2}(i,1),2), handles.newvertices(handles.class{2}(i,1),3), 1, 0, 0);  
                  end
                  for i = 1:size(handles.class{3},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{3}(i,1),1), handles.newvertices(handles.class{3}(i,1),2), handles.newvertices(handles.class{3}(i,1),3), 0, 1, 0);  
                  end
                  for i = 1:size(handles.class{4},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{4}(i,1),1), handles.newvertices(handles.class{4}(i,1),2), handles.newvertices(handles.class{4}(i,1),3), 1, 1, 0);  
                  end 
                  for i = 1:size(handles.class{5},1)
                      fprintf(fileID,'v %f %f %f %d %d %d\n',handles.newvertices(handles.class{5}(i,1),1), handles.newvertices(handles.class{5}(i,1),2), handles.newvertices(handles.class{5}(i,1),3), 1, 1, 0);  
                  end                     
               fclose(fileID);  
               
               newsubmesh1 = strcat(B{1},'_1.txt');
              fileID = fopen(newsubmesh1, 'w');
                for i = 1:size(handles.class{1},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{1}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh2 = strcat(B{1},'_2.txt');
              fileID = fopen(newsubmesh2, 'w');
                for i = 1:size(handles.class{2},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{2}(i,1),1));  
                end
              fclose(fileID);
              
              newsubmesh3 = strcat(B{1},'_3.txt');
              fileID = fopen(newsubmesh3, 'w');
                for i = 1:size(handles.class{3},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{3}(i,1),1));  
                end
              fclose(fileID);               

              newsubmesh4 = strcat(B{1},'_4.txt');
              fileID = fopen(newsubmesh4, 'w');
                for i = 1:size(handles.class{4},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{4}(i,1),1));  
                end
              fclose(fileID); 

              newsubmesh5 = strcat(B{1},'_5.txt');
              fileID = fopen(newsubmesh5, 'w');
                for i = 1:size(handles.class{5},1)
                    fprintf(fileID,'%d \t',red_features(handles.class{5}(i,1),1));  
                end
              fclose(fileID);               
           end
         end

    end
    
    message.String = ' '; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Function for changing the size of the points using slider %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = slider_change_size(varargin)

    if  handles.ok == 0
        
        gui_pop()

    else    

        cla(ax)
        cla(ax1)
        cla(ax2)

        message.String = 'In processing !!!'; 

        % Callback for the edit box and slider.
        [h,Slider] = varargin{[1,3]};  % Get calling handle and structure.
        SL = get(Slider.sl,{'min','value','max'});  % Get the slider's info.
        newvalue = get(Slider.ed(2),{'String'}); 
        E = floor(str2double(get(h,'string')));  % Numerical edit string.

        switch h  % Who called?
            case Slider.ed(1)
                if E <= SL{2}
                    set(Slider.sl,'min',E)  % E is less than current value.
                elseif E < SL{3}
                    set(Slider.sl,'val',E,'min',E) % E is less than max value.
                    set(Slider.ed(2),'string',E) % Set the current display.
                else
                    set(h,'string',SL{1}) % Reset the value.
                end
            case Slider.ed(2)
                if E >= SL{1} && E <= SL{3}
                    set(Slider.sl,'value',E)  % E falls within range of slider.
                else
                    set(h,'string',SL{2}) % User tried to set slider out of range. 
                end
            case Slider.ed(3)
                if E >= SL{2}
                    set(Slider.sl,'max',E)  % E is less than current value.
                elseif E > SL{1}
                    set(Slider.sl,'val',E,'max',E) % E is less than max value.
                    set(Slider.ed(2),'string',E) % Set the current display.
                else
                    set(h,'string',SL{3}) % Reset the value.
                end      
            case Slider.sl
                set(Slider.ed(2),'string',SL{2}) % Set edit to current slider.
            otherwise
                % Do nothing
        end

        % Change the size of the displayed points
        axes(ax)
        scatter3(handles.Original(:,1), handles.Original(:,2), handles.Original(:,3), str2double(newvalue), 'b', 'filled')
        rotate3d on
        axis tight

        axes(ax1)
        scatter3(handles.blue(:,1), handles.blue(:,2), handles.blue(:,3), str2double(newvalue), 'b', 'filled')
        hold on
        scatter3(handles.red(:,1), handles.red(:,2), handles.red(:,3), str2double(newvalue), 'r', 'filled')
        rotate3d on
        axis tight
        
        
        axes(ax2)
        if handles.kBest==1
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
        end
        if handles.kBest==2
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
        end
        if handles.kBest==3
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on
        end
        if handles.kBest==4
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{4},1), handles.newvertices(handles.class{4},2), handles.newvertices(handles.class{4},3), str2double(newvalue), 'y', 'filled')
            hold on
        end
        if handles.kBest==5
            scatter3(handles.newvertices(handles.class{1},1), handles.newvertices(handles.class{1},2), handles.newvertices(handles.class{1},3), str2double(newvalue), 'b', 'filled')
            hold on 
            scatter3(handles.newvertices(handles.class{2},1), handles.newvertices(handles.class{2},2), handles.newvertices(handles.class{2},3), str2double(newvalue), 'r', 'filled')
            hold on
            scatter3(handles.newvertices(handles.class{3},1), handles.newvertices(handles.class{3},2), handles.newvertices(handles.class{3},3), str2double(newvalue), 'g', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{4},1), handles.newvertices(handles.class{4},2), handles.newvertices(handles.class{4},3), str2double(newvalue), 'y', 'filled')
            hold on            
            scatter3(handles.newvertices(handles.class{5},1), handles.newvertices(handles.class{5},2), handles.newvertices(handles.class{5},3), str2double(newvalue), 'c', 'filled')
            hold on
        end

        message.String = ' ';
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function for the pop up window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = gui_pop()
    
    msgbox('You need to load a .PLY file first','Warning Message')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
