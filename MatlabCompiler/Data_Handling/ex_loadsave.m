function ex_loadsave
% This example shows how to work with the "load/save" functions
% on data files in deployed mode. There are three source data files
% in this example.
%    user_data.mat 
%    userdata/extra_data.mat 
%    ../externdata/extern_data.mat
%
% Compile this example with mcc command: 
%     mcc -m ex_loadsave.m -a 'user_data.mat'
%         -a './userdata/extra_data.mat'
%         -a '../externdata/extern_data.mat'
%
% All '-a' data files will be included in the CTF archive and expanded under ctfroot.
% User data files are unchanged by mcc. There is no excryption of user included data files.
%
% Any '-a' file located in or under the main '-m' file's directory at compile time,
% will be located at the same path relative to ctfroot at run time. In this example,
% the run time location of 'user_data.mat' will be '$ctfroot\user_data.mat',
% and that of 'extra_data.mat' will be '$ctfroot\userdata\extra_data.mat'.
%
% Any '-a' file NOT located in or under the main '-m' file's directory at compile time,
% will have its absolute path from the disk drive or filesystem root appended to
% ctfroot to give its run time location. In this example, since the compile time
% location of 'extern_data.mat' is 'c:\$matlabroot\examples\externdata\extern_data.mat',
% its run time location will be '$ctfroot\$matlabroot\examples\externdata\extern_data.mat'.
%
% The output data file, relative to the application's run time current working directory, is:
%   ./output/saved_data.mat
% When writing data files to local disk, do not save any files under ctfroot,
% as ctfroot may be deleted and/or refreshed when the application is next started.
%
%==== load data file =============================
if isdeployed
    % In deployed mode, all files in or under the main file's directory
    % may be loaded by full path, or by path relative to ctfroot.
    % LOADFILENAME1=which(fullfile(ctfroot,'user_data.mat'));    
    % LOADFILENAME2=which(fullfile(ctfroot,'userdata','extra_data.mat'));
    LOADFILENAME1=which(fullfile('user_data.mat'));
    LOADFILENAME2=which(fullfile('userdata','extra_data.mat'));
    % For external data file, full path will be added into ctf, you don't
    % need specify the full path to find the file.
    LOADFILENAME3=which(fullfile('extern_data.mat'));
else
    %running the code in MATLAB
    LOADFILENAME1=fullfile(matlabroot,'extern','examples','compiler','Data_Handling','user_data.mat');
    LOADFILENAME2=fullfile(matlabroot,'extern','examples','compiler','Data_Handling','userdata','extra_data.mat');
    LOADFILENAME3=fullfile(matlabroot,'extern','examples','compiler','externdata','extern_data.mat');
end

% Load the data file from current working directory
disp(['Load A from : ',LOADFILENAME1]);
load(LOADFILENAME1,'data1');
disp('A= ');
disp(data1);

% Load the data file from sub directory
disp(['Load B from : ',LOADFILENAME2]);
load(LOADFILENAME2,'data2');
disp('B= ');
disp(data2);

% Load extern data outside of current working directory
disp(['Load extern data from : ',LOADFILENAME3]);
load(LOADFILENAME3);
disp('ext_data= ');
disp(ext_data);

%==== multiple the data matrix by 2 ==============
result = data1*data2;
disp('A * B = ');
disp(result);

%==== save  the new data to a new file ===========
SAVEPATH=strcat(pwd,filesep,'output');
if ( ~isdir(SAVEPATH))
    mkdir(SAVEPATH);
end
SAVEFILENAME=strcat(SAVEPATH,filesep,'saved_data.mat');
disp(['Save the A * B result to : ',SAVEFILENAME]);
save(SAVEFILENAME, 'result');




