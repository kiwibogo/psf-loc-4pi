


datapath = 'D:\Sheng\data\11-06-2021 DM scan\dm676\';
filelist = dir([datapath,'*3phase*.mat']);

fn = length(filelist);
for ii = 1:3
    
    load([datapath,filelist(ii).name])
    imsz = size(qd1);
    for ss = 1:4
        eval(['qd',num2str(ss),' = reshape(qd',num2str(ss),',imsz(1),imsz(2),[]);'])
    end
    save(filelist(ii).name,'qd1','qd2','qd3','qd4')
end
