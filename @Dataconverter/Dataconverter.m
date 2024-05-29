classdef Dataconverter < handle
    properties
        Maindir;
        Folderlist;
        Copydir;
    end
    
    methods
        function obj = Dataconverter()
            obj.gui;
        end
      
        
        function dcimg2mat_4pi(obj,type)
            dirN = length(obj.Folderlist);
            hp = waitbar(0,'converting data');
            for ii=1:dirN
                if iscell(obj.Folderlist)
                    datapath = obj.Folderlist{ii};
                else
                    datapath = fullfile(obj.Maindir,obj.Folderlist(ii).name);
                end
                datafiles = dir([datapath,'\', '*.dcimg']);
                %centers =[388, 573, 1400, 1586];
                %centers = [282,490,1673,1881];
                %centers = [368, 535, 1724, 1891];
                %centers = [391,558,1699,1868];
                %centers = [391,557,1702,1869];
                %centers = [344,511,1748,1916];
                %centers = [362,529,1730,1898];
                %centers = [362,529,1726,1893];
                %centers = [123,356,1583,1815];
                %centers = [174,431,1620,1877];
                %centers = [175,432,1619,1876];
                centers = dlmread('C:\Users\Ries Lab\Documents\sheng-gitlab\4Pi-analysis-SL\quad_center.txt','\t');

                switch type
                    case 'concatenate'
                        dcimg2mat(datapath,datafiles,centers,1)
                    case 'none'
                        dcimg2mat(datapath,datafiles,centers,0)
                end
                waitbar(ii/dirN,hp);
            end
            
            close(hp)
        end
        
        function dcimg2mat_palm2d(obj,type)
            dirN = length(obj.Folderlist);
            hp = waitbar(0,'converting data');
            for ii=1:dirN
                if iscell(obj.Folderlist)
                    datapath = obj.Folderlist{ii};
                else
                    datapath = fullfile(obj.Maindir,obj.Folderlist(ii).name);
                end
                datafiles = dir([datapath,'\', '*.dcimg']);
                switch type
                    case 'concatenate'
                        dcimg2mat(datapath,datafiles,[],1)
                    case 'none'
                        dcimg2mat(datapath,datafiles,[],0)
                end
                waitbar(ii/dirN,hp);
            end
            close(hp)
        end
        
        function selectfiles(obj,keyword)
            if isempty(obj.Maindir) %|| (obj.Maindir==0)
                obj.Maindir = uigetdir('');
            else
                obj.Maindir = uigetdir(obj.Maindir);
            end
            list1 = dir([obj.Maindir,'\',keyword,'*']);
            obj.Folderlist = list1(cell2mat({list1.isdir}));
            obj.rmconvertfolder();
        end
        
        function rmconvertfolder(obj)
            list1 = obj.Folderlist;
            fn = length(list1);
            ind = 1;
            for ii = 1:fn
                if ind>length(list1)
                    break;
                end
                files = dir([list1(ind).folder,'\', list1(ind).name,'\*.dcimg']);
                if isempty(files)
                    list1(ind) = [];
                else
                    ind = ind+1;
                end
            end
            obj.Folderlist = list1;
        end
        
        function copyfiles(obj)
            dirN = length(obj.Folderlist);
            hp = waitbar(0,'copy data');
            for ii = 1:dirN
                if iscell(obj.Folderlist)
                    datapath = obj.Folderlist{ii};
                else
                    datapath = fullfile(obj.Maindir,obj.Folderlist(ii).name);
                end
                
                if ~exist(obj.Copydir,'dir')
                    mkdir(obj.Copydir)
                end
                
                cd([datapath,'\'])
                if ~isempty(dir('*.mat'))
                    copyfile('*.mat',obj.Copydir);
                end
                waitbar(ii/dirN,hp);
            end
            close(hp)
        end
        
        function deletedcimg(obj)
            dirN = length(obj.Folderlist);
            hp = waitbar(0,'delete dcimg files');
            for ii = 1:dirN
                if iscell(obj.Folderlist)
                    datapath = obj.Folderlist{ii};
                else
                    datapath = fullfile(obj.Maindir,obj.Folderlist(ii).name);
                end
                
                delete([datapath,'\*.dcimg'])
                waitbar(ii/dirN,hp);
            end
            close(hp)
        end
        
    end
end