function gui(obj)
xsz = 530;
ysz = 280;
xst = 200;
yst = 300;

guiFig = figure('Units','pixels','Position',[xst yst xsz ysz],...
    'MenuBar','none','ToolBar','figure','Visible','on',...
    'NumberTitle','off','UserData',0,'Tag',...
    'Dataconverter.gui','HandleVisibility','off','name','convert to mat files',...
    'CloseRequestFcn',@figclose);

defaultBackground = get(0,'defaultUicontrolBackgroundColor');
set(guiFig,'Color',defaultBackground);
handles.output = guiFig;
guidata(guiFig,handles);

handles.list_file = uicontrol('Parent',guiFig,'Style','listbox','String','','Value',1,'Units','Normalized','Position',[.03 .5 .7 .4]);
handles.edit_key = uicontrol('Parent',guiFig,'Style','edit','String','cell','Fontsize',10,'Units','Normalized','Position',[.75 .65 .22 .1]);
handles.button_selectf = uicontrol('Parent',guiFig, 'Style', 'pushbutton', 'String','select files','Enable','on','Units','Normalized','Position', [.77 .5 .15 .1],'Fontsize',10,'Callback',@selectf);
handles.popup_type = uicontrol('Parent',guiFig, 'Style', 'popupmenu', 'String',{'concatenate','none'},'Enable','on','Value',2,'Units','Normalized','Position', [.03 .35 .15 .06],'Fontsize',10,'Callback',@updatefolder);
handles.popup_datatype = uicontrol('Parent',guiFig, 'Style', 'popupmenu', 'String',{'4pi','palm 2d'},'Enable','on','Value',1,'Units','Normalized','Position', [.2 .35 .15 .06],'Fontsize',10);

handles.button_conv = uicontrol('Parent',guiFig, 'Style', 'pushbutton', 'String','convert to mat','Enable','on','Units','Normalized','Position', [.42 .32 .2 .1],'Fontsize',10,'Callback',@conv2mat);
handles.edit_cpdir = uicontrol('Parent',guiFig,'Style','edit','String','','Fontsize',10,'Units','Normalized','Position',[.03 .15 .7 .1]);
handles.button_copyf = uicontrol('Parent',guiFig, 'Style', 'pushbutton', 'String','copy files','Enable','on','Units','Normalized','Position', [.77 .03 .15 .1],'Fontsize',10,'Callback',@copyf);
handles.button_setf = uicontrol('Parent',guiFig, 'Style', 'pushbutton', 'String','set copy dir','Enable','on','Units','Normalized','Position', [.77 .15 .15 .1],'Fontsize',10,'Callback',@setdir);

cfile = uicontextmenu(guiFig);
handles.list_file.UIContextMenu = cfile;
uimenu('Parent',cfile,'Label','delete dcimg','Callback',@deleteraw);

    function figclose(~,~)
        delete(guiFig)
    end

    function selectf(~,~)
        keyword = get(handles.edit_key,'String');
        obj.selectfiles(keyword)
        N1 = length(obj.Folderlist);
        folderlist = cell(1,N1);
        for ii = 1:N1
            folderlist{ii} = fullfile(obj.Maindir,obj.Folderlist(ii).name);
        end
        set(handles.list_file,'String',folderlist);  
        %[~,foldername] = fileparts(obj.Maindir);
        obj.Copydir = obj.Maindir;
        set(handles.edit_cpdir,'String',obj.Copydir);
    end

    function conv2mat(~,~)
        type = get(handles.popup_type,'String');
        ind = get(handles.popup_type,'Value');
        dtype = get(handles.popup_datatype,'String');
        ind1 = get(handles.popup_datatype,'Value');
        switch dtype{ind1}
            case '4pi'
                obj.dcimg2mat_4pi(type{ind})
            case 'palm 2d'
                obj.dcimg2mat_palm2d(type{ind})
        end
    end
    
    function updatefolder(h,~)
        N1 = length(obj.Folderlist);
        folderlist = cell(1,N1);
        for ii = 1:N1
            folderlist{ii} = fullfile(obj.Maindir,obj.Folderlist(ii).name);
        end
        set(handles.list_file,'String',folderlist);  
    end
    
    function setdir(~,~)
        obj.Copydir = uigetdir('Y:\');
        set(handles.edit_cpdir,'String',obj.Copydir);
    end

    function copyf(~,~)
        obj.Copydir = get(handles.edit_cpdir,'String');
        obj.copyfiles();
    end

    function deleteraw(~,~)
        obj.deletedcimg();
    end

end