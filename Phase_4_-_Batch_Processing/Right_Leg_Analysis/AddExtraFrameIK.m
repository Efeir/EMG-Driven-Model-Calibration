function AddExtraFrameIK
% This function adds an extra time frame to the end of a marker motion file
% because opensim will often ignore the final time frame during inverse
% kinematics

files = dir('*.trc');
nfiles = length(files);

for i = 1:nfiles

infile = files(i).name
    
fid = fopen(infile);
header1 = fgetl(fid);
header2 = fgetl(fid);
header3 = fgetl(fid);
header4 = fgetl(fid);
header5 = fgetl(fid);

% Read marker data as one long column
data = fscanf(fid,'%f');
fclose(fid);

% Reshape marker data based on number of rows and
% columns of marker data specified in the input data
% file, with two extra columns for frame and time
info = sscanf(header3,'%f %f %d %d');
nrows = info(3,1);
ncols = info(4,1)*3;
data = reshape(data,ncols+2,nrows)';

% Separate frame and time vectors from marker data
frame_time_data = data(:,1:2);
marker_data = data(:,3:ncols+2);

step = mean(diff(frame_time_data));
frames_time_data = [frame_time_data; frame_time_data(end,:)+step];
marker_data = [marker_data; marker_data(end,:)];

time = frames_time_data(:,2);

MarkerNames = strrep(header4,'Patient 4:','');
MarkerNames = textscan(MarkerNames,'%s');
MarkerNames = MarkerNames{1};
MarkerNames(1:2) = [];

outputTRCFile(infile, marker_data , time,MarkerNames)

end

function outputTRCFile(outFile, MarkerData,TimeVector,markerNames)

nMarkers = length(markerNames);

Header = HeaderTRC;

rate = 1/mean(diff(TimeVector));
AllMarkersTRC = MarkerData;

MarkerLabels = markerNames;
Header = strrep(Header,'DATARATE_INS',sprintf('%.12f',rate));
Header = strrep(Header,'CAMRATE_INS',sprintf('%.12f',rate));
Header = strrep(Header,'ORIGRATE_INS',sprintf('%.12f',rate));
Header = strrep(Header,'NFRAMES_INS',sprintf('%d',length(TimeVector)));
Header = strrep(Header,'NMARKERS_INS',sprintf('%d',nMarkers));
Header = strrep(Header,'UNITS_INS','mm');
Header = strrep(Header,'FILENAME',outFile);

% Create Marker Labels String
RemoveXYZLabel = 0;
if length(MarkerLabels)==(nMarkers*3)
    MarkerLabels=MarkerLabels([1:3:end]);
    RemoveXYZLabel = 1;
end
MarkerLabelString ='';
XYZString ='';
for m=1:nMarkers
    if length(MarkerLabels)<m
        label = ['Marker' num2str(m)];
    else
        label = MarkerLabels{m};
        if RemoveXYZLabel
            label = strrep(label,'_x','');
        end
    end
    XYZString = [XYZString sprintf('X%d\tY%d\tZ%d\t',m,m,m)];
    MarkerLabelString = [MarkerLabelString sprintf('%s\t\t\t',label)];
end
Header = strrep(Header,'FIRSTMARKERNAME_INS',MarkerLabelString);
Header = strrep(Header,'XYZ_INS',XYZString);

%% Create Data string
DataString = '';
for f=1:length(TimeVector)
    DataString = [DataString sprintf('%d\t%.12f\t',f,TimeVector(f))];
    for m=1:nMarkers
        DataString = [DataString sprintf('%.12f\t%.12f\t%.12f\t',AllMarkersTRC(f,(m-1)*3+1),AllMarkersTRC(f,(m-1)*3+2),AllMarkersTRC(f,(m-1)*3+3))];
    end
    DataString = [DataString sprintf('\n')];
end
Header = strrep(Header,'DATASTART_INS',DataString);

fid = fopen(outFile,'w');
fprintf(fid,'%s',Header);
fclose(fid);

function header = HeaderTRC

header = ['PathFileType	4	(X/Y/Z)	FILENAME.trc' sprintf('\n') ...
    'DataRate	CameraRate	NumFrames	NumMarkers	Units	OrigDataRate	OrigDataStartFrame	OrigNumFrames' sprintf('\n') ...
    'DATARATE_INS	CAMRATE_INS	NFRAMES_INS	NMARKERS_INS	UNITS_INS	ORIGRATE_INS	1	       NFRAMES_INS' sprintf('\n') ...
    'Frame#	Time	FIRSTMARKERNAME_INS' sprintf('\n') ...
    sprintf('\t')  sprintf('\t')	'XYZ_INS' sprintf('\n') ...
    '' sprintf('\n') ...
    'DATASTART_INS' sprintf('\n')];