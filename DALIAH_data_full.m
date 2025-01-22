    %% Load data
% Change path and save old path
oldFolder = cd('M:\data_cancer\DALIAH\Trine_DALIAH_5y');

% Load data
dataTable = readtable('j_ottesen_180822 - Tobias');
dataTableExtra = readtable('j_ottesen_new_variables_180822 - Tobias');
dataTableExtraDates = readtable('edatemed_randomized_LOCKED(220509)');

% Change path back
cd(oldFolder);

% Dimensions
[M,N] = size(dataTable);
P = M;

% Remove patients not used in large data set from pre-baseline
temp = ismember(dataTableExtra.id,dataTable.id);
dataTableExtra(~temp,:) = [];
temp = ismember(dataTableExtraDates.id,dataTable.id);
dataTableExtraDates(~temp,:) = [];

%% Data processing
% Change weird observations
dataTable.b21fmd(2) = 4;
dataTable.b25ohaem(108) = 8.2;
dataTable.b21bdate(194) = datetime('10-jun-2017'); % midpoint between the two other available dates
dataTable.b21bdate(200) = datetime('18-oct-2018'); % midpoint between the two other available dates

% Number of measurements (including pre-baseline)
nMeasurements = 26;

% Cohort data
mut = dataTableExtra.b01_drivermutation;
endVisit = dataTableExtra.edateend_visit;
endTreatVisit = dataTableExtraDates.edatemed_randomized; % Updated dates from Trine
priorHU = dataTableExtra.b00hu;
priorHUDur = dataTableExtra.b00hu_sum_days;
priorHUTotDose = dataTableExtra.b00hu_sum_dose;
assignedTreatmentType = dataTableExtra.tx_allocation;
hospID = dataTable.id;
diagnosis = dataTable.b01bwho;
age = dataTable.b01bage;
gender = dataTable.b01bsex;
initspleen = dataTable.b01bstr;
initspleen2 = dataTable.b01bsul1;
pruritus = dataTable.b01akloe;

% Throw away some extra updated dates

% Streamline mutations
mut(strcmp(mut,'CALR type 2')) = {'CALR'};
mut(strcmp(mut,'CALR, type 1')) = {'CALR'};
mut(strcmp(mut,'CALR, variant 34 bp del')) = {'CALR'};
mut(strcmp(mut,'MPL, W515K')) = {'MPL'};
mut(strcmp(mut,'MPL, unspecified' )) = {'MPL'};

% Convert priorHU to logical
temp = false(size(priorHU));
temp(strcmp(priorHU,'Yes')) = 1;
priorHU = temp;

% Count mutations
[dlistMut,iaMut,icMut] = unique(mut);
dcountsMut = accumarray(icMut,1);
dlistCountsMut = {dlistMut, dcountsMut};

% Binary variable for male
temp = false(size(gender));
temp(strcmp(gender,'Male')) = 1;
malep = temp;

% Convert pruritus to binary
temp = false(size(pruritus));
temp(strcmp(pruritus,'Yes')) = 1;
pruritus = temp;

% Storage
patientID = NaN(nMeasurements*M,1);
baselineID = false(nMeasurements*M,1);
dates = NaT(nMeasurements*M,1);
relTime = NaN(nMeasurements*M,1);
HGB = NaN(nMeasurements*M,1);
WBC = NaN(nMeasurements*M,1);
neu = NaN(nMeasurements*M,1);
HCT = NaN(nMeasurements*M,1);
TRC = NaN(nMeasurements*M,1);
jak = NaN(nMeasurements*M,1);
phlebot = NaN(nMeasurements*M,1);
RBCTrans = NaN(nMeasurements*M,1);
treatType = strings(nMeasurements*M,1);
treatUnits = strings(nMeasurements*M,1);
treatDose = NaN(nMeasurements*M,1);
treat2Type = strings(nMeasurements*M,1);
treat2Units = strings(nMeasurements*M,1);
treat2Dose = NaN(nMeasurements*M,1);
IFNsys1 = NaN(nMeasurements*M,1);
IFNsys2 = NaN(nMeasurements*M,1);
IFNintron1 = NaN(nMeasurements*M,1);
IFNintron2 = NaN(nMeasurements*M,1);
HU1 = NaN(nMeasurements*M,1);
HU2 = NaN(nMeasurements*M,1);
HU3 = NaN(nMeasurements*M,1);

% Assign visit numbers
visitNumber = repmat((0:25)',[M,1]);

% Assign IDs in loop
for i=1:M
    % Assign ID
    patientID((i-1)*nMeasurements+1:i*nMeasurements) = i;
end

% Assign baseline ID
baselineID(2:nMeasurements:end) = 1;

% Load pre-baseline data
tempstr = "b0" + 0;
dates(1:nMeasurements:end) = dataTableExtra.(tempstr + "bdate");
relTime(1:nMeasurements:end) = days(dataTableExtra.(tempstr + "bdate")...
                              - dataTable.b01bdate);
HGB(1:nMeasurements:end) = dataTableExtra.(tempstr + "ohaem");
WBC(1:nMeasurements:end) = dataTableExtra.(tempstr + "oleu");
neu(1:nMeasurements:end) = dataTableExtra.(tempstr + "oneu");
HCT(1:nMeasurements:end) = dataTableExtra.(tempstr + "oef");

% Load baseline data
tempstr = "b0" + 1;
dates(2:nMeasurements:end) = dataTable.(tempstr + "bdate");
relTime(2:nMeasurements:end) = days(dataTable.(tempstr + "bdate")...
                              - dataTable.b01bdate);
HGB(2:nMeasurements:end) = dataTable.(tempstr + "ohaem");
WBC(2:nMeasurements:end) = dataTable.(tempstr + "oleu");
neu(2:nMeasurements:end) = dataTable.(tempstr + "oneu");
TRC(2:nMeasurements:end) = dataTable.(tempstr + "otrc");
HCT(2:nMeasurements:end) = dataTable.(tempstr + "oef");
phlebot(2:nMeasurements:end) = dataTable.(tempstr + "avene");
RBCTrans(2:nMeasurements:end) = dataTable.(tempstr + "asagm");
tempTreatType = string(dataTable.(tempstr + "fffm"));
tempTreatUnits = string(dataTable.(tempstr + "fmd"));

% Since we have no other information, the initial dose is taken from the
% current dose at the second measurement
tempTreatDose = dataTable.b02bfmbd;

% Extra treatment also seems to be typed in wrongly. I assume that a number
% in dose corresponds to Hydrea.
tempTreat2Dose = dataTable.(tempstr + "fsbd");
tempTreat2Type = strings(size(tempTreat2Dose));
tempTreat2Units = strings(size(tempTreat2Dose));
tempTreat2Type(~isnan(tempTreat2Dose)) = "Hydrea";
tempTreat2Units(~isnan(tempTreat2Dose)) = "mg/day";

% Divide doses to get per day and adjust units
tempTreatDose(strcmp(tempTreatUnits,"mg/2 days")) = tempTreatDose(strcmp(tempTreatUnits,"mg/2 days"))/2;
tempTreatDose(strcmp(tempTreatUnits,"mg/2. day")) = tempTreatDose(strcmp(tempTreatUnits,"mg/2. day"))/2;
tempTreatDose(strcmp(tempTreatUnits,"mg/7 days")) = tempTreatDose(strcmp(tempTreatUnits,"mg/7 days"))/7;
tempTreatDose(strcmp(tempTreatUnits,"mg/week")) = tempTreatDose(strcmp(tempTreatUnits,"mg/week"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/7 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/7 days"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/week")) = tempTreatDose(strcmp(tempTreatUnits,"ug/week"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/10 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/10 days"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/14 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/14 days"))/14;
tempTreatDose(strcmp(tempTreatUnits,"ug/21 days")) =  tempTreatDose(strcmp(tempTreatUnits,"ug/21 days"))/21;
tempTreatDose(strcmp(tempTreatUnits,"ug/28 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/28 days"))/28;
tempTreatDose(strcmp(tempTreatUnits,"ug/35 days")) =  tempTreatDose(strcmp(tempTreatUnits,"ug/35 days"))/35;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/2 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/2 days"))/2;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/2. day")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/2. day"))/2;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/7 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/7 days"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/week")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/week"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/7 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/7 days"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/week")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/week"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/10 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/10 days"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/14 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/14 days"))/14;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/21 days")) =  tempTreat2Dose(strcmp(tempTreat2Units,"ug/21 days"))/21;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/28 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/28 days"))/28;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/35 days")) =  tempTreat2Dose(strcmp(tempTreat2Units,"ug/35 days"))/35;
tempTreatUnits(strcmp(tempTreatUnits,"mg/2 days")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"mg/2. day")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"mg/7 days")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"mg/week")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/7 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/week")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/10 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/14 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/21 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/28 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/35 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/2 days")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/2. day")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/7 days")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/week")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/7 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/week")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/10 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/14 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/21 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/28 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/35 days")) = "ug/day";

% Assign temporary variables to permanent variables
treatType(2:nMeasurements:end) = tempTreatType;
treatUnits(2:nMeasurements:end) = tempTreatUnits;
treatDose(2:nMeasurements:end) = tempTreatDose;
treat2Type(2:nMeasurements:end) = tempTreat2Type;
treat2Units(2:nMeasurements:end) = tempTreat2Units;
treat2Dose(2:nMeasurements:end) = tempTreat2Dose;

% Load post-baseline data (except last visit), treatments as temporary
for i=2:nMeasurements-2
    % Different strings when below/above 10
    if i<9
        tempstr = "b0" + i;
        tempstr2 = "b0" + (i+1);
    elseif i==9
        tempstr = "b0" + i;
        tempstr2 = "b" + (i+1);    
    elseif i>9
        tempstr = "b" + i;
        tempstr2 = "b" + (i+1);
    end
    dates(i+1:nMeasurements:end) = dataTable.(tempstr + "bdate");
    relTime(i+1:nMeasurements:end) = days(dataTable.(tempstr + "bdate")...
                                  - dataTable.b01bdate);
    HGB(i+1:nMeasurements:end) = dataTable.(tempstr + "ohaem");
    WBC(i+1:nMeasurements:end) = dataTable.(tempstr + "oleu");
    neu(i+1:nMeasurements:end) = dataTable.(tempstr + "oneu");
    TRC(i+1:nMeasurements:end) = dataTable.(tempstr + "otrc");
    HCT(i+1:nMeasurements:end) = dataTable.(tempstr + "oef");
    phlebot(i+1:nMeasurements:end) = dataTable.(tempstr + "avene");
    RBCTrans(i+1:nMeasurements:end) = dataTable.(tempstr + "asagm");
    tempTreatType = string(dataTable.(tempstr2 + "bfm"));
    tempTreatUnits = string(dataTable.(tempstr2 + "bfmde"));
    tempTreat2Type = string(dataTable.(tempstr2 + "bsbeh"));
    tempTreat2Units = string(dataTable.(tempstr2 + "bsbde"));

    % Extract backward doses from next measurement
    tempTreatDose = dataTable.(tempstr2 + "bfmbd");
    tempTreat2Dose = dataTable.(tempstr2 + "bsbd");

    % Streamline dose units
    tempTreatUnits(strcmp(tempTreatUnits,"1")) = "mg/day";
    tempTreatUnits(strcmp(tempTreatUnits,"2")) = "mg/2 days";
    tempTreatUnits(strcmp(tempTreatUnits,"3")) = "mg/7 days";
    tempTreatUnits(strcmp(tempTreatUnits,"4")) = "ug/7 days";
    tempTreatUnits(strcmp(tempTreatUnits,"5")) = "ug/10 days";
    tempTreatUnits(strcmp(tempTreatUnits,"6")) = "ug/14 days";
    tempTreatUnits(strcmp(tempTreatUnits,"7")) = "ug/21 days";
    tempTreatUnits(strcmp(tempTreatUnits,"8")) = "ug/28 days";
    tempTreatUnits(strcmp(tempTreatUnits,"9")) = "ug/35 days";
    tempTreat2Units(strcmp(tempTreat2Units,"1")) = "mg/day";
    tempTreat2Units(strcmp(tempTreat2Units,"2")) = "mg/2 days";
    tempTreat2Units(strcmp(tempTreat2Units,"3")) = "mg/7 days";
    tempTreat2Units(strcmp(tempTreat2Units,"4")) = "ug/7 days";
    tempTreat2Units(strcmp(tempTreat2Units,"5")) = "ug/10 days";
    tempTreat2Units(strcmp(tempTreat2Units,"6")) = "ug/14 days";
    tempTreat2Units(strcmp(tempTreat2Units,"7")) = "ug/21 days";
    tempTreat2Units(strcmp(tempTreat2Units,"8")) = "ug/28 days";
    tempTreat2Units(strcmp(tempTreat2Units,"9")) = "ug/35 days";
    
    % Divide doses to get per day and adjust units
    tempTreatDose(strcmp(tempTreatUnits,"mg/2 days")) = tempTreatDose(strcmp(tempTreatUnits,"mg/2 days"))/2;
    tempTreatDose(strcmp(tempTreatUnits,"mg/2. day")) = tempTreatDose(strcmp(tempTreatUnits,"mg/2. day"))/2;
    tempTreatDose(strcmp(tempTreatUnits,"mg/7 days")) = tempTreatDose(strcmp(tempTreatUnits,"mg/7 days"))/7;
    tempTreatDose(strcmp(tempTreatUnits,"mg/week")) = tempTreatDose(strcmp(tempTreatUnits,"mg/week"))/7;
    tempTreatDose(strcmp(tempTreatUnits,"ug/7 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/7 days"))/7;
    tempTreatDose(strcmp(tempTreatUnits,"ug/week")) = tempTreatDose(strcmp(tempTreatUnits,"ug/week"))/7;
    tempTreatDose(strcmp(tempTreatUnits,"ug/10 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/10 days"))/7;
    tempTreatDose(strcmp(tempTreatUnits,"ug/14 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/14 days"))/14;
    tempTreatDose(strcmp(tempTreatUnits,"ug/21 days")) =  tempTreatDose(strcmp(tempTreatUnits,"ug/21 days"))/21;
    tempTreatDose(strcmp(tempTreatUnits,"ug/28 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/28 days"))/28;
    tempTreatDose(strcmp(tempTreatUnits,"ug/35 days")) =  tempTreatDose(strcmp(tempTreatUnits,"ug/35 days"))/35;
    tempTreat2Dose(strcmp(tempTreat2Units,"mg/2 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/2 days"))/2;
    tempTreat2Dose(strcmp(tempTreat2Units,"mg/2. day")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/2. day"))/2;
    tempTreat2Dose(strcmp(tempTreat2Units,"mg/7 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/7 days"))/7;
    tempTreat2Dose(strcmp(tempTreat2Units,"mg/week")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/week"))/7;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/7 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/7 days"))/7;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/week")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/week"))/7;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/10 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/10 days"))/7;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/14 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/14 days"))/14;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/21 days")) =  tempTreat2Dose(strcmp(tempTreat2Units,"ug/21 days"))/21;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/28 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/28 days"))/28;
    tempTreat2Dose(strcmp(tempTreat2Units,"ug/35 days")) =  tempTreat2Dose(strcmp(tempTreat2Units,"ug/35 days"))/35;
    tempTreatUnits(strcmp(tempTreatUnits,"mg/2 days")) = "mg/day";
    tempTreatUnits(strcmp(tempTreatUnits,"mg/2. day")) = "mg/day";
    tempTreatUnits(strcmp(tempTreatUnits,"mg/7 days")) = "mg/day";
    tempTreatUnits(strcmp(tempTreatUnits,"mg/week")) = "mg/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/7 days")) = "ug/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/week")) = "ug/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/10 days")) = "ug/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/14 days")) = "ug/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/21 days")) = "ug/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/28 days")) = "ug/day";
    tempTreatUnits(strcmp(tempTreatUnits,"ug/35 days")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"mg/2 days")) = "mg/day";
    tempTreat2Units(strcmp(tempTreat2Units,"mg/2. day")) = "mg/day";
    tempTreat2Units(strcmp(tempTreat2Units,"mg/7 days")) = "mg/day";
    tempTreat2Units(strcmp(tempTreat2Units,"mg/week")) = "mg/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/7 days")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/week")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/10 days")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/14 days")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/21 days")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/28 days")) = "ug/day";
    tempTreat2Units(strcmp(tempTreat2Units,"ug/35 days")) = "ug/day";

    % Assign doses etc. from temporary variables to real variables
    treatType(i+1:nMeasurements:end) = tempTreatType;
    treatUnits(i+1:nMeasurements:end) = tempTreatUnits;
    treatDose(i+1:nMeasurements:end) = tempTreatDose;
    treat2Type(i+1:nMeasurements:end) = tempTreat2Type;
    treat2Units(i+1:nMeasurements:end) = tempTreat2Units;
    treat2Dose(i+1:nMeasurements:end) = tempTreat2Dose;
end

% Load last visit
tempstr = "b" + 25;
dates(26:nMeasurements:end) = dataTable.(tempstr + "bdate");
relTime(26:nMeasurements:end) = days(dataTable.(tempstr + "bdate")...
                              - dataTable.b01bdate);
HGB(26:nMeasurements:end) = dataTable.(tempstr + "ohaem");
WBC(26:nMeasurements:end) = dataTable.(tempstr + "oleu");
neu(26:nMeasurements:end) = dataTable.(tempstr + "oneu");
TRC(26:nMeasurements:end) = dataTable.(tempstr + "otrc");
HCT(26:nMeasurements:end) = dataTable.(tempstr + "oef");
phlebot(26:nMeasurements:end) = dataTable.(tempstr + "avene");
RBCTrans(26:nMeasurements:end) = dataTable.(tempstr + "asagm");
tempTreatType = string(dataTable.(tempstr + "bfm"));
tempTreat2Type = string(dataTable.(tempstr + "bsbeh"));
tempTreatUnits = string(dataTable.(tempstr + "bfmde"));
tempTreat2Units = string(dataTable.(tempstr + "bsbde"));

% Since we have no other information, the dose is assumed to continue
tempTreatDose = dataTable.(tempstr + "bfmbd");
tempTreat2Dose = dataTable.(tempstr + "bsbd");

% Streamline dose units
tempTreatUnits(strcmp(tempTreatUnits,"1")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"2")) = "mg/2 days";
tempTreatUnits(strcmp(tempTreatUnits,"3")) = "mg/7 days";
tempTreatUnits(strcmp(tempTreatUnits,"4")) = "ug/7 days";
tempTreatUnits(strcmp(tempTreatUnits,"5")) = "ug/10 days";
tempTreatUnits(strcmp(tempTreatUnits,"6")) = "ug/14 days";
tempTreatUnits(strcmp(tempTreatUnits,"7")) = "ug/21 days";
tempTreatUnits(strcmp(tempTreatUnits,"8")) = "ug/28 days";
tempTreatUnits(strcmp(tempTreatUnits,"9")) = "ug/35 days";
tempTreat2Units(strcmp(tempTreat2Units,"1")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"2")) = "mg/2 days";
tempTreat2Units(strcmp(tempTreat2Units,"3")) = "mg/7 days";
tempTreat2Units(strcmp(tempTreat2Units,"4")) = "ug/7 days";
tempTreat2Units(strcmp(tempTreat2Units,"5")) = "ug/10 days";
tempTreat2Units(strcmp(tempTreat2Units,"6")) = "ug/14 days";
tempTreat2Units(strcmp(tempTreat2Units,"7")) = "ug/21 days";
tempTreat2Units(strcmp(tempTreat2Units,"8")) = "ug/28 days";
tempTreat2Units(strcmp(tempTreat2Units,"9")) = "ug/35 days";

% Divide doses to get per day and adjust units
tempTreatDose(strcmp(tempTreatUnits,"mg/2 days")) = tempTreatDose(strcmp(tempTreatUnits,"mg/2 days"))/2;
tempTreatDose(strcmp(tempTreatUnits,"mg/2. day")) = tempTreatDose(strcmp(tempTreatUnits,"mg/2. day"))/2;
tempTreatDose(strcmp(tempTreatUnits,"mg/7 days")) = tempTreatDose(strcmp(tempTreatUnits,"mg/7 days"))/7;
tempTreatDose(strcmp(tempTreatUnits,"mg/week")) = tempTreatDose(strcmp(tempTreatUnits,"mg/week"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/7 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/7 days"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/week")) = tempTreatDose(strcmp(tempTreatUnits,"ug/week"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/10 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/10 days"))/7;
tempTreatDose(strcmp(tempTreatUnits,"ug/14 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/14 days"))/14;
tempTreatDose(strcmp(tempTreatUnits,"ug/21 days")) =  tempTreatDose(strcmp(tempTreatUnits,"ug/21 days"))/21;
tempTreatDose(strcmp(tempTreatUnits,"ug/28 days")) = tempTreatDose(strcmp(tempTreatUnits,"ug/28 days"))/28;
tempTreatDose(strcmp(tempTreatUnits,"ug/35 days")) =  tempTreatDose(strcmp(tempTreatUnits,"ug/35 days"))/35;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/2 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/2 days"))/2;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/2. day")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/2. day"))/2;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/7 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/7 days"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"mg/week")) = tempTreat2Dose(strcmp(tempTreat2Units,"mg/week"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/7 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/7 days"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/week")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/week"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/10 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/10 days"))/7;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/14 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/14 days"))/14;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/21 days")) =  tempTreat2Dose(strcmp(tempTreat2Units,"ug/21 days"))/21;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/28 days")) = tempTreat2Dose(strcmp(tempTreat2Units,"ug/28 days"))/28;
tempTreat2Dose(strcmp(tempTreat2Units,"ug/35 days")) =  tempTreat2Dose(strcmp(tempTreat2Units,"ug/35 days"))/35;
tempTreatUnits(strcmp(tempTreatUnits,"mg/2 days")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"mg/2. day")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"mg/7 days")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"mg/week")) = "mg/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/7 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/week")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/10 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/14 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/21 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/28 days")) = "ug/day";
tempTreatUnits(strcmp(tempTreatUnits,"ug/35 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/2 days")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/2. day")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/7 days")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"mg/week")) = "mg/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/7 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/week")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/10 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/14 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/21 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/28 days")) = "ug/day";
tempTreat2Units(strcmp(tempTreat2Units,"ug/35 days")) = "ug/day";

% Assign temporary variables to permanent variables
treatType(26:nMeasurements:end) = tempTreatType;
treatUnits(26:nMeasurements:end) = tempTreatUnits;
treatDose(26:nMeasurements:end) = tempTreatDose;
treat2Type(26:nMeasurements:end) = tempTreat2Type;
treat2Units(26:nMeasurements:end) = tempTreat2Units;
treat2Dose(26:nMeasurements:end) = tempTreat2Dose;

% Load drop-out times as months
DOTime = days(dataTable.edateend- dataTable.b01bdate)/30.4;

% Load allele burden measurements
jak(1+1:nMeasurements:end) = dataTable.b01ojak_validated;
jak(5+1:nMeasurements:end) = dataTable.b05ojak_validated;
jak(7+1:nMeasurements:end) = dataTable.b07ojak_validated;
jak(9+1:nMeasurements:end) = dataTable.b09ojak_validated;
jak(11+1:nMeasurements:end) = dataTable.b11ojak_validated;
jak(13+1:nMeasurements:end) = dataTable.b13ojak_validated;
jak(17+1:nMeasurements:end) = dataTable.b17ojak_validated;
jak(21+1:nMeasurements:end) = dataTable.b21ojak_validated;
jak(25+1:nMeasurements:end) = dataTable.b25ojak_validated;

% ID for different treatments
IDsys = strcmp(treatType,'Pegasys');
ID2sys = strcmp(treat2Type,'Yes, Pegasys');
IDintron = strcmp(treatType,'PegIntron');
ID2intron = strcmp(treat2Type,'Yes, PegIntron');
IDHU = strcmp(treatType,'Hydrea');
ID2HU = strcmp(treat2Type,'Yes, Hydrea');
ID3HU = strcmp(treat2Type,'Hydrea');

% ID for different diagnoses
PVIDp = strcmp(diagnosis,'Polycythemia vera');
ETIDp = strcmp(diagnosis,'Essential thrombocythemia');
MFID1p = strcmp(diagnosis,'Prefibrotic myelofibrosis');
MFID2p = strcmp(diagnosis,'Primary myelofibrosis');
MFIDp = logical(MFID1p+MFID2p);

% Storage
PVID = false(M*nMeasurements,1);
ETID = false(M*nMeasurements,1);
MFID = false(M*nMeasurements,1);

% Assign observations to diagnosis IDs
for i=1:P
    % Assign to right group
    if PVIDp(i)
        PVID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
    elseif ETIDp(i)
        ETID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
    elseif MFIDp(i)
        MFID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
    end
end

% Put doses in correct treatment variables
IFNsys1(IDsys) = treatDose(IDsys);
IFNsys2(ID2sys) = treat2Dose(ID2sys);
IFNintron1(IDintron) = treatDose(IDintron);
IFNintron2(ID2intron) = treat2Dose(ID2intron);
HU1(IDHU) = treatDose(IDHU);
HU2(ID2HU) = treat2Dose(ID2HU);
HU3(ID3HU) = treat2Dose(ID3HU);

% Adjust IFN values that lack a unit
tempID = logical(ismissing(treatUnits) + (treatUnits == ""));
tempID2 = logical(ismissing(treat2Units) + (treat2Units == ""));
IFNsys1(tempID) = IFNsys1(tempID)/7;
IFNsys2(tempID2) = IFNsys2(tempID2)/7;
IFNintron1(tempID) = IFNintron1(tempID)/7;
IFNintron2(tempID2) = IFNintron2(tempID2)/7;

% Collect different treatments
IFNsys = sum([IFNsys1,IFNsys2],2,'omitnan');
IFNintron = sum([IFNintron1,IFNintron2],2,'omitnan');
HU = sum([HU1,HU2,HU3],2,'omitnan');
IFN = IFNsys+IFNintron;

% Convert haematocrit and jak to decimal
HCT = HCT/100;
jak = jak/100;

% Convert time to months
relTime = relTime/30.4;

% Choose time limits for plots (in months)
minRelTime = -20;
maxRelTime = 65;

% Storage
IFNIDp = false(P,1);
sysIDp = false(P,1);
intronIDp = false(P,1);
sysIntronIDp = false(P,1);
HUIDp = false(P,1);
MixIDp = false(P,1);
IFNID = false(M*nMeasurements,1);
sysID = false(M*nMeasurements,1);
intronID = false(M*nMeasurements,1);
sysIntronID = false(M*nMeasurements,1);
HUID = false(M*nMeasurements,1);
MixID = false(M*nMeasurements,1);

% Assign patients to treatment groups
for i=1:P
    % Define ID for the patient
    pID = i;
    IDvec = patientID==pID;

    % Temporary logicals
    tempHUID = sum(HU(IDvec))>0;
    tempIFNID1 = sum(IFNsys(IDvec))>0;
    tempIFNID2 = sum(IFNintron(IDvec))>0;

    % Assign to right group
    if tempHUID && ~tempIFNID1 && ~tempIFNID2
        HUIDp(i) = true;
        HUID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
    elseif ~tempHUID && ( tempIFNID1 || tempIFNID2 )
        IFNIDp(i) = true;
        IFNID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
        if tempIFNID1 && ~tempIFNID2
            sysIDp(i) = true;
            sysID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
        elseif ~tempIFNID1 && tempIFNID2
            intronIDp(i) = true;
            intronID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
        else
            sysIntronIDp(i) = true;
            sysIntronID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
        end
    else
        MixIDp(i) = true;
        MixID(1+(i-1)*nMeasurements:i*nMeasurements) = true;
    end
end

% Dates registered
timeStudyVisitCat = [-1; 0; 0.5; 1; (2:2:12)'; (15:3:63)'];

% Storage
timeStudyVisit = nan(M*nMeasurements,1);

% Assign register dates in loop
for i=1:M*nMeasurements
    % Find minimum difference
    temp = abs(relTime(i)-timeStudyVisitCat);
    [minTemp,ID] = min(temp);
    
   if ~isnan(minTemp)
       % Assign date
       timeStudyVisit(i) = timeStudyVisitCat(ID);
   end
end

% % Check for wrong catgories
% catDiff = diff(timeStudyVisit);
% 
% % Find IDs for wrong categories
% ID = catDiff==0;
% 
% % Adjust ID's
% for i=1:M-1
%     if ID(i)
%         ID2 = find(timeStudyVisit(i+1)==timeStudyVisitCat);
%         timeStudyVisit(i+1) = timeStudyVisitCat(ID2+1);
%     end
% end


% Relative allele burdens
% Storage
relJak = zeros(size(jak));
basJak = zeros(P,1);

% Calculate relative allele burden
for i=1:P
    % Find patient
    pID = i;
    ID = patientID==pID;
    
    % Find baseline
    if ~isempty(jak(logical(ID.*baselineID)))
        basJak(i) = jak(logical(ID.*baselineID));
    else
        basJak(i) = NaN;
    end

    % Calculate relative allele burden
    relJak(ID) = jak(ID)/basJak(i);
end

% Number of observations and eligible patients for fitting
% Storage
numObs = zeros(P,1);
maxJak = NaN(P,1);

% Calculate number of observations in loop
for i=1:P
    pID = i;
    IDvec = patientID==pID;
    numObs(i) = sum(~isnan(jak(IDvec)));
end

% Storage
eliID = false(M,1);

% Make eligible ID, 3 or more observations
for i=1:M*nMeasurements
    % Extract patientID
    j = patientID(i);
    if j>0
        % Check if 3 or more observations
        eliID(i) = numObs(j)>=3;
    end
end

% Patient-specific eligible ID
pEliID = numObs>=3;

%% Numbers for IFN patients with prior HU treatment
% Extract ID for patients fitted
fitIDp = logical(pEliID.*IFNIDp);

%% Check number of phlebotomies
% Storage
numPhlebot = zeros(P,1);
for i=1:P
    pID = i;
    tempID = patientID == pID;
    numPhlebot = sum(phlebot(tempID),'omitnan');
end

%% Number checks
% Unique diagnoses with counts
[dlist,ia,ic] = unique(diagnosis);
dcounts = accumarray(ic,1);
dlistCounts = {dlist, dcounts};

% Choose observation times to include for drop-outs
inclTime = [0;4;8;12;18;24;36;48;60];

% Storage
nParti = zeros(size(inclTime));
nPartiPV = zeros(size(inclTime));
nPartiET = zeros(size(inclTime));
nPartiMF = zeros(size(inclTime));

% Calculate participants at included times
for i=1:length(inclTime)
    nParti(i) = sum(DOTime-(inclTime(i))>-1);
    nPartiPV(i) = sum(DOTime(PVIDp)-(inclTime(i))>-1);
    nPartiET(i) = sum(DOTime(ETIDp)-(inclTime(i))>-1);
    nPartiMF(i) = sum(DOTime(MFIDp)-(inclTime(i))>-1);
end

% Plot
figure();
hold on;
grid on;
plot(inclTime,nParti/max(nParti),'.--k','linewidth',1.5,'markersize',20);
plot(inclTime,nPartiPV/max(nParti),'.--r','linewidth',1.5,'markersize',20);
plot(inclTime,nPartiET/max(nParti),'.--g','linewidth',1.5,'markersize',20);
plot(inclTime,nPartiMF/max(nParti),'.--b','linewidth',1.5,'markersize',20);
title('Fraction of Patients over Time','interpreter','latex','fontsize',16);
xlabel('Time/months','interpreter','latex','fontsize',12);
ylabel('Fraction','interpreter','latex','fontsize',12);
xlim([0,max(inclTime)]);
ylim([0,1]);
legend({'Total','PV','ET','MF'},'fontsize',12,'location','best',...
        'interpreter','latex');

%% Plot allele burden for all patients in the same plot
% Figure
figure();

% Loop
for i=1:P
    % Define ID for the patient
    pID = i;
    IDvec = patientID==pID;
    
    % Find the baseline date and hospital ID
    basDate = dates(logical(IDvec.*baselineID));
    hospIDp = hospID(pID);
    
    % Vectors for plotting
    time = relTime(IDvec);
    jakVec = jak(IDvec);
    
    % Vector for treatment start
    baselinex = days(linspace(basDate,basDate,1001)-basDate)/30.4;
    baseliney = linspace(0,100,1001);
    
    % Plot
    hold on;
    % plot(days(preDates(pID)-basDate)/30.4,prejak(pID),'.r','MarkerSize',16);
    plot(time(~isnan(jakVec)),jakVec(~isnan(jakVec))*100,'.-','MarkerSize',16,'linewidth',1.5);
    % plot(baselinex,baseliney,'-k','linewidth',1);
end

% Title, labels etc.
title('\textit{JAK2} VAF for All Patients - DALIAH',...
          'interpreter','latex','fontsize',16);
xlabel('Time/months','interpreter','latex','fontsize',12);
ylabel('\textit{JAK2} VAF/\%','interpreter','latex','fontsize',12);
xlim([0,maxRelTime]);
ylim([0,100]);

% Save figure
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('jakAllD1')),'jpg');
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('jakAllD1')),'epsc');

% Figure
figure();

% Loop
for i=1:P
    % Define ID for the patient
    pID = i;
    IDvec = patientID==pID;
    
    % Find the baseline date and hospital ID
    basDate = dates(logical(IDvec.*baselineID));
    hospIDp = hospID(pID);
    
    % Vectors for plotting
    time = relTime(IDvec);
    relJakVec = relJak(IDvec);
    
    % Vector for treatment start
    baselinex = days(linspace(basDate,basDate,1001)-basDate)/30.4;
    baseliney = linspace(0,100,1001);
    
    % Plot
    hold on;
    % plot(days(preDates(pID)-basDate)/30.4,prejak(pID),'.r','MarkerSize',16);
    plot(time(~isnan(relJakVec)),relJakVec(~isnan(relJakVec)),'.-','MarkerSize',16,'linewidth',1.5);
    % plot(baselinex,baseliney,'-k','linewidth',1);
end

% Title, labels etc.
title('Relative \textit{JAK2} VAF for All Patients - DALIAH',...
          'interpreter','latex','fontsize',16);
xlabel('Time/months','interpreter','latex','fontsize',12);
ylabel('$\frac{JAK2}{JAK2_0}$/1','interpreter','latex','fontsize',12);
xlim([0,maxRelTime]);
ylim([0,1.5]);

% Save figure
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('relJakAllD1')),'jpg');
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('relJakAllD1')),'epsc');

%% Plot allele burden for all patients in the same plot
% Figure
figure();

% Loop
for i=1:P
    % Define ID for the patient
    pID = i;
    IDvec = patientID==pID;
    
    % Find the baseline date and hospital ID
    basDate = dates(logical(IDvec.*baselineID));
    hospIDp = hospID(pID);
    
    % Vectors for plotting
    time = relTime(IDvec);
    jakVec = jak(IDvec);
    
    % Vector for treatment start
    baselinex = days(linspace(basDate,basDate,1001)-basDate)/30.4;
    baseliney = linspace(0,100,1001);
    
    % Plot
    hold on;
    % plot(days(preDates(pID)-basDate)/30.4,prejak(pID),'.r','MarkerSize',16);
    plot(time(~isnan(jakVec)),jakVec(~isnan(jakVec))*100,'.-','MarkerSize',16,'linewidth',1.5);
    % plot(baselinex,baseliney,'-k','linewidth',1);
end

% Title, labels etc.
title('Allele Burden for All Patients - DALIAH',...
          'interpreter','latex','fontsize',16);
xlabel('Time/months','interpreter','latex','fontsize',12);
ylabel('Allele burden','interpreter','latex','fontsize',12);
xlim([0,maxRelTime]);
ylim([0,100]);

% Save figure
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('jakAllD')),'jpg');
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('jakAllD')),'epsc');

% Figure
figure();

% Loop
for i=1:P
    % Define ID for the patient
    pID = i;
    IDvec = patientID==pID;
    
    % Find the baseline date and hospital ID
    basDate = dates(logical(IDvec.*baselineID));
    hospIDp = hospID(pID);
    
    % Vectors for plotting
    time = relTime(IDvec);
    relJakVec = relJak(IDvec);
    
    % Vector for treatment start
    baselinex = days(linspace(basDate,basDate,1001)-basDate)/30.4;
    baseliney = linspace(0,100,1001);
    
    % Plot
    hold on;
    % plot(days(preDates(pID)-basDate)/30.4,prejak(pID),'.r','MarkerSize',16);
    plot(time(~isnan(relJakVec)),relJakVec(~isnan(relJakVec)),'.-','MarkerSize',16,'linewidth',1.5);
    % plot(baselinex,baseliney,'-k','linewidth',1);
end

% Title, labels etc.
title('Allele Burden for All Patients - DALIAH',...
          'interpreter','latex','fontsize',16);
xlabel('Time/months','interpreter','latex','fontsize',12);
ylabel('Allele burden','interpreter','latex','fontsize',12);
xlim([0,maxRelTime]);
ylim([0,1.5]);

% Save figure
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('relJakAllD')),'jpg');
% saveas(gcf,fullfile('H:\Kode\Billeder og video\', sprintf('relJakAllD')),'epsc');

%% Allele burden medians
% Eligibility of patients for JAK2 on both patient and observation level,
% also observations to the end
eliJakp = false(P,1);
eliJak = false(M,1);
for p=1:P
    % Patient ID
    pID = p;
    IDvec = patientID==pID;
    
    % Assign
    if sum(jak(IDvec)>0)>0 % && max(timeStudyVisit(IDvec)>21)
        eliJakp(p) = 1;
        eliJak(IDvec) = 1;
    end
end


% Choose observation times to include
inclTime = [0;4;8;12;18;24;36;48;60];

% Storage
jakMedsIFN = zeros(size(inclTime));
jakMedsIFNCI = zeros(size(jakMedsIFN,1),2);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) median(x,'omitnan');

% Calculate medians at different times and confidence intervals
for i=1:length(jakMedsIFN)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*IFNID);
    jakMedsIFN(i) = median(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMedsIFNCI(i,:) = tempCI';
end

% Plot
figure();
hold on;
plot(inclTime,jakMedsIFN*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMedsIFNCI(:,1)'*100 fliplr(jakMedsIFNCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMedsIFNCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMedsIFNCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Median Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Median allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample median','Bootstrap CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
jakMedsHU = zeros(size(inclTime));
jakMedsHUCI = zeros(size(jakMedsHU,1),2);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) median(x,'omitnan');

% Calculate medians at different times and confidence intervals
for i=1:length(jakMedsHU)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*HUID);
    jakMedsHU(i) = median(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMedsHUCI(i,:) = tempCI';
end

% Plot
figure();
hold on;
plot(inclTime,jakMedsHU*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMedsHUCI(:,1)'*100 fliplr(jakMedsHUCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMedsHUCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMedsHUCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Median Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Median allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample median','Bootstrap CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
jakMedsMix = zeros(size(inclTime));
jakMedsMixCI = zeros(size(jakMedsMix,1),2);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) median(x,'omitnan');

% Calculate medians at different times and confidence intervals
for i=1:length(jakMedsMix)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*MixID);
    jakMedsMix(i) = median(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMedsMixCI(i,:) = tempCI';
end

% Plot
figure();
hold on;
plot(inclTime,jakMedsMix*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMedsMixCI(:,1)'*100 fliplr(jakMedsMixCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMedsMixCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMedsMixCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Median Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Median allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample median','Bootstrap CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
jakMedsPV = zeros(size(inclTime));
jakMedsPVCI = zeros(size(jakMedsPV,1),2);
jakMedsET = zeros(size(inclTime));
jakMedsETCI = zeros(size(jakMedsET,1),2);
jakMedsMF = zeros(size(inclTime));
jakMedsMFCI = zeros(size(jakMedsMF,1),2);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) median(x,'omitnan');

% Calculate medians at different times and confidence intervals
for i=1:length(jakMedsPV)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*PVID.*IFNID);
    jakMedsPV(i) = median(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMedsPVCI(i,:) = tempCI';
    % jakMEDSPVsd()
    ID = logical(temp.*ETID.*IFNID);
    jakMedsET(i) = median(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMedsETCI(i,:) = tempCI';
    ID = logical(temp.*MFID.*IFNID);
    jakMedsMF(i) = median(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMedsMFCI(i,:) = tempCI';
end

% Plot
figure();
hold on;
plot(inclTime,jakMedsPV*100,'.-r','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMedsPVCI(:,1)'*100 fliplr(jakMedsPVCI(:,2)'*100)], ...
      'r', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMedsPVCI(:,1)*100,'-r','linewidth',1.5);
plot(inclTime,jakMedsPVCI(:,2)*100,'-r','linewidth',1.5);
plot(inclTime,jakMedsET*100,'.-g','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMedsETCI(:,1)'*100 fliplr(jakMedsETCI(:,2)'*100)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMedsETCI(:,1)*100,'-g','linewidth',1.5);
plot(inclTime,jakMedsETCI(:,2)*100,'-g','linewidth',1.5);
plot(inclTime,jakMedsMF*100,'.-b','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMedsMFCI(:,1)'*100 fliplr(jakMedsMFCI(:,2)'*100)], ...
      'b', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMedsMFCI(:,1)*100,'-b','linewidth',1.5);
plot(inclTime,jakMedsMFCI(:,2)*100,'-b','linewidth',1.5);
ylim([0,100]);
title('Median Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Median allele burden/\%','fontsize',12,'interpreter','latex');
legend({'PV - Sample median','PV - Bootstrap CI','','',...
        'ET - Sample median','ET - Bootstrap CI','','',...
        'MF - Sample median','MF - Bootstrap CI','',''},...
        'fontsize',12,'location','best','interpreter','latex');

%% Boxplots of relative allele burdens and fit to medians - "fit"
% Choose observation times to include
inclTime = [0;4;8;12;18;24;36;48;60];

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) median(x,'omitnan');

% Storage
relJakIFNMeds = zeros(size(inclTime));
relJakIFNMedsCI = zeros(size(relJakIFNMeds,1),2);

% Calculate medians at different times
for i=1:length(relJakIFNMeds)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*eliID.*IFNID);
    relJakIFNMeds(i) = median(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakIFNMedsCI(i,:) = tempCI';
end

% Storage
relJakHUMeds = zeros(size(inclTime));
relJakHUMedsCI = zeros(size(relJakHUMeds,1),2);

% Calculate medians at different times
for i=1:length(relJakHUMeds)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*eliID.*HUID);
    relJakHUMeds(i) = median(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakHUMedsCI(i,:) = tempCI';
end

% Storage
relJakMixMeds = zeros(size(inclTime));
relJakMixMedsCI = zeros(size(relJakMixMeds,1),2);

% Calculate medians at different times
for i=1:length(relJakMixMeds)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*eliID.*MixID);
    relJakMixMeds(i) = median(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMixMedsCI(i,:) = tempCI';
end

% Locations
% xloc = unique(timeStudyVisit);
xloc = inclTime;

% Calculate exponential fit using "fit"
fo = fitoptions('Method','NonlinearLeastSquares','Lower',0, ...
                'Upper',Inf,'StartPoint',1);
ft = fittype('exp(-a*x)','dependent',{'y'},'independent',{'x'},...
             'coefficients',{'a'},'options',fo);
[expfit,expgof] = fit(inclTime,relJakIFNMeds,ft);

% Calculate biexponential fit
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0], ...
                'Upper',[Inf,Inf,Inf],'StartPoint',[1,1,1]);
ft = fittype(['B*((beta^2+c^2+0.49)/c^2*exp(-beta^2*x)-(beta^2+0.49)' ...
              '/c^2*exp(-(beta^2+c^2)*x))'],'dependent',{'y'},...
             'independent',{'x'},'coefficients',{'B','beta','c'},...
             'options',fo);
[bifit,bigof] = fit(inclTime,relJakIFNMeds,ft);

% Select only relevant months for boxplot
ID = ismember(timeStudyVisit,inclTime);
timeStudyVisitBox = timeStudyVisit(ID);
relJakBox = relJak(ID);

% Boxplot
figure();
hold on;
boxplot(relJakBox,timeStudyVisitBox,'Positions',xloc);

% Plot fitted curve
plot(bifit,inclTime,relJakIFNMeds);

% Boxplot
figure();
hold on;
boxplot(relJakBox,timeStudyVisitBox,'Positions',xloc);

% Plot fitted curve
plot(expfit,inclTime,relJakIFNMeds);

% Plot confidence and prediction intervals
CI = predint(expfit,inclTime,0.95,'functional','on');
PI = predint(expfit,inclTime,0.95,'observation','on');
plot(inclTime,CI(:,1),'g');
plot(inclTime,PI(:,1),'c');
plot(inclTime,CI(:,2),'g');
plot(inclTime,PI(:,2),'c');

% Labels etc.
title('Fit using MATLAB''s "fit"','fontsize',16,'interpreter','latex');
xlabel('Time since baseline/months','fontsize',12,'interpreter','latex');
ylabel('$\frac{JAK2}{JAK2_0}$','fontsize',16,'interpreter','latex');
legend({'Medians','Exponential fit','Confidence interval',...
        'Prediction interval'},'location','best','fontsize',12);

% Storage
relJakMedsPV = zeros(size(inclTime));
relJakMedsPVCI = zeros(size(relJakMedsPV,1),2);
relJakMedsET = zeros(size(inclTime));
relJakMedsETCI = zeros(size(relJakMedsET,1),2);
relJakMedsMF = zeros(size(inclTime));
relJakMedsMFCI = zeros(size(relJakMedsMF,1),2);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) median(x,'omitnan');

% Calculate medians at different times and confidence intervals
for i=1:length(relJakMedsPV)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*PVID.*IFNID);
    relJakMedsPV(i) = median(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMedsPVCI(i,:) = tempCI';
    ID = logical(temp.*ETID.*IFNID);
    relJakMedsET(i) = median(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMedsETCI(i,:) = tempCI';
    ID = logical(temp.*MFID.*IFNID);
    relJakMedsMF(i) = median(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMedsMFCI(i,:) = tempCI';
end

% Plot
figure();
hold on;
plot(inclTime,relJakMedsPV,'.-r','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMedsPVCI(:,1)' fliplr(relJakMedsPVCI(:,2)')], ...
      'r', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMedsPVCI(:,1),'-r','linewidth',1.5);
plot(inclTime,relJakMedsPVCI(:,2),'-r','linewidth',1.5);
plot(inclTime,relJakMedsET,'.-g','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMedsETCI(:,1)' fliplr(relJakMedsETCI(:,2)')], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMedsETCI(:,1),'-g','linewidth',1.5);
plot(inclTime,relJakMedsETCI(:,2),'-g','linewidth',1.5);
plot(inclTime,relJakMedsMF,'.-b','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMedsMFCI(:,1)' fliplr(relJakMedsMFCI(:,2)')], ...
      'b', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMedsMFCI(:,1),'-b','linewidth',1.5);
plot(inclTime,relJakMedsMFCI(:,2),'-b','linewidth',1.5);
ylim([0,1]);
title('Median Relative Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Median allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample median','Bootstrap CI'},'fontsize',12,'location','best',...
        'interpreter','latex');
legend({'PV - Sample relative median','PV - Bootstrap CI','','',...
        'ET - Sample relative median','ET - Bootstrap CI','','',...
        'MF - Sample relative median','MF - Bootstrap CI','',''},...
        'fontsize',12,'location','best','interpreter','latex');

%% Allele burden means
% Eligibility of patients for JAK2 on both patient and observation level,
% also observations to the end
eliJakp = false(P,1);
eliJak = false(M,1);
for p=1:P
    % Patient ID
    pID = p;
    IDvec = patientID==pID;
    
    % Assign
    if sum(jak(IDvec)>0)>0 % && max(timeStudyVisit(IDvec)>21)
        eliJakp(p) = 1;
        eliJak(IDvec) = 1;
    end
end


% Choose observation times to include
inclTime = [0;4;8;12;18;24;36;48;60];

% Storage
jakMeansIFN = zeros(size(inclTime));
jakMeansIFNCI = zeros(size(jakMeansIFN,1),2);
jakMeansIFNCIstd = zeros(size(jakMeansIFN,1),2);
nIFN = zeros(size(jakMeansIFN));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(jakMeansIFN)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*IFNID);
    jakMeansIFN(i) = mean(jak(ID),'omitnan');
    temp2 = jak(ID);
    n = length(temp2(~isnan(temp2)));
    nIFN(i) = n;
    tempCI = bootci(nboot,{meanfun,jak(ID)},'Alpha',0.05);
    jakMeansIFNCI(i,:) = tempCI';
    tempCI = [jakMeansIFN(i)-tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n);...
              jakMeansIFN(i)+tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n)];
    jakMeansIFNCIstd(i,:) = tempCI';
end


% Plot
figure();
hold on;
plot(inclTime,jakMeansIFN*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansIFNCI(:,1)'*100 fliplr(jakMeansIFNCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
patch([inclTime' fliplr(inclTime')], [jakMeansIFNCIstd(:,1)'*100 fliplr(jakMeansIFNCIstd(:,2)'*100)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansIFNCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMeansIFNCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
jakMeansHU = zeros(size(inclTime));
jakMeansHUCI = zeros(size(jakMeansHU,1),2);
nHU = zeros(size(jakMeansHU));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(jakMeansHU)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*HUID);
    jakMeansHU(i) = mean(jak(ID),'omitnan');
    temp2 = jak(ID);
    n = length(temp2(~isnan(temp2)));
    nHU(i) = n;
%     tempCI = [jakMeansHU(i)-tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n);...
%               jakMeansHU(i)+tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n)];
%     jakMeansHUCI(i,:) = tempCI';
    tempCI = bootci(nboot,{meanfun,jak(ID)},'Alpha',0.05);
    jakMeansHUCI(i,:) = tempCI';
end


% Plot
figure();
hold on;
plot(inclTime,jakMeansHU*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansHUCI(:,1)'*100 fliplr(jakMeansHUCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansHUCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMeansHUCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
jakMeansMix = zeros(size(inclTime));
jakMeansMixCI = zeros(size(jakMeansMix,1),2);
nMix = zeros(size(jakMeansMix));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(jakMeansMix)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*MixID);
    jakMeansMix(i) = mean(jak(ID),'omitnan');
    temp2 = jak(ID);
    n = length(temp2(~isnan(temp2)));
    nMix(i) = n;
%     tempCI = [jakMeansMix(i)-tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n);...
%               jakMeansMix(i)+tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n)];
%     jakMeansMixCI(i,:) = tempCI';
    tempCI = bootci(nboot,{meanfun,jak(ID)},'Alpha',0.05);
    jakMeansMixCI(i,:) = tempCI';
end


% Plot
figure();
hold on;
plot(inclTime,jakMeansMix*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansMixCI(:,1)'*100 fliplr(jakMeansMixCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansMixCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMeansMixCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
jakMeansPV = zeros(size(inclTime));
jakMeansPVCI = zeros(size(jakMeansPV,1),2);
jakMeansPVsd = zeros(size(jakMeansPV,1),1);
jakMeansET = zeros(size(inclTime));
jakMeansETCI = zeros(size(jakMeansET,1),2);
jakMeansETsd = zeros(size(jakMeansPV,1),1);
jakMeansMF = zeros(size(inclTime));
jakMeansMFCI = zeros(size(jakMeansMF,1),2);
jakMeansMFsd = zeros(size(jakMeansPV,1),1);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) mean(x,'omitnan');

% Calculate means at different times and confidence intervals
for i=1:length(jakMeansPV)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*PVID.*IFNID);
    jakMeansPV(i) = mean(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMeansPVCI(i,:) = tempCI';
    jakMeansPVsd(i) = std(jak(ID),'omitnan');
    ID = logical(temp.*ETID.*IFNID);
    jakMeansET(i) = mean(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMeansETCI(i,:) = tempCI';
    jakMeansETsd(i) = std(jak(ID),'omitnan');
    ID = logical(temp.*MFID.*IFNID);
    jakMeansMF(i) = mean(jak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,jak(ID)},'Alpha',0.05);
    jakMeansMFCI(i,:) = tempCI';
    jakMeansMFsd(i) = std(jak(ID),'omitnan');
end

% Plot
figure();
hold on;
plot(inclTime,jakMeansPV*100,'.-r','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansPVCI(:,1)'*100 fliplr(jakMeansPVCI(:,2)'*100)], ...
      'r', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansPVCI(:,1)*100,'-r','linewidth',1.5);
plot(inclTime,jakMeansPVCI(:,2)*100,'-r','linewidth',1.5);
plot(inclTime,jakMeansET*100,'.-g','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansETCI(:,1)'*100 fliplr(jakMeansETCI(:,2)'*100)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansETCI(:,1)*100,'-g','linewidth',1.5);
plot(inclTime,jakMeansETCI(:,2)*100,'-g','linewidth',1.5);
plot(inclTime,jakMeansMF*100,'.-b','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansMFCI(:,1)'*100 fliplr(jakMeansMFCI(:,2)'*100)], ...
      'b', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansMFCI(:,1)*100,'-b','linewidth',1.5);
plot(inclTime,jakMeansMFCI(:,2)*100,'-b','linewidth',1.5);
ylim([0,100]);
title('Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'PV - Sample mean','PV - Bootstrap CI','','',...
        'ET - Sample mean','ET - Bootstrap CI','','',...
        'MF - Sample mean','MF - Bootstrap CI','',''},...
        'fontsize',12,'location','best','interpreter','latex');

figure();
hold on;
errorbar(inclTime,jakMeansPV*100,jakMeansPVsd*100,'.-r','linewidth',1.5,'markersize',20);
errorbar(inclTime,jakMeansET*100,jakMeansETsd*100,'.-g','linewidth',1.5,'markersize',20);
errorbar(inclTime,jakMeansMF*100,jakMeansMFsd*100,'.-b','linewidth',1.5,'markersize',20);
ylim([0,100]);
title('Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'PV - Sample mean $\pm$ standard deviation'
        'ET - Sample mean $\pm$ standard deviation'
        'MF - Sample mean $\pm$ standard deviation',},...
        'fontsize',12,'location','best','interpreter','latex');

%% Allele burden means for 64 patients eligible for data fitting and treated with IFN
% Eligibility of patients for JAK2 on both patient and observation level,
% also observations to the end
eliJakp = false(P,1);
eliJak = false(M,1);
for p=1:P
    % Patient ID
    pID = p;
    IDvec = patientID==pID;
    
    % Assign
    if sum(jak(IDvec)>0)>0 % && max(timeStudyVisit(IDvec)>21)
        eliJakp(p) = 1;
        eliJak(IDvec) = 1;
    end
end


% Choose observation times to include
inclTime = [0;4;8;12;18;24;36;48;60];

% Storage
jakMeansIFNEli = zeros(size(inclTime));
jakMeansIFNEliCI = zeros(size(jakMeansIFNEli,1),2);
jakMeansIFNEliCIstd = zeros(size(jakMeansIFNEli,1),2);
nIFN = zeros(size(jakMeansIFNEli));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(jakMeansIFNEli)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*IFNID.*eliID);
    jakMeansIFNEli(i) = mean(jak(ID),'omitnan');
    temp2 = jak(ID);
    n = length(temp2(~isnan(temp2)));
    nIFN(i) = n;
    tempCI = bootci(nboot,{meanfun,jak(ID)},'Alpha',0.05);
    jakMeansIFNEliCI(i,:) = tempCI';
    tempCI = [jakMeansIFNEli(i)-tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n);...
              jakMeansIFNEli(i)+tinv(0.975,n-1)*std(jak(ID),'omitnan')/sqrt(n)];
    jakMeansIFNEliCIstd(i,:) = tempCI';
    min(jak(ID))
    max(jak(ID))
end


% Plot
figure();
hold on;
plot(inclTime,jakMeansIFNEli*100,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [jakMeansIFNEliCI(:,1)'*100 fliplr(jakMeansIFNEliCI(:,2)'*100)], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
patch([inclTime' fliplr(inclTime')], [jakMeansIFNEliCIstd(:,1)'*100 fliplr(jakMeansIFNEliCIstd(:,2)'*100)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,jakMeansIFNEliCI(:,1)*100,'-k','linewidth',1.5);
plot(inclTime,jakMeansIFNEliCI(:,2)*100,'-k','linewidth',1.5);
ylim([0,100]);
title('Temporal Evolution of Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

%% Allele burden mean relative
% Eligibility of patients for JAK2 on both patient and observation level,
% also observations to the end
eliJakp = false(P,1);
eliJak = false(M,1);
for p=1:P
    % Patient ID
    pID = p;
    IDvec = patientID==pID;
    
    % Assign
    if sum(relJak(IDvec)>0)>0 % && max(timeStudyVisit(IDvec)>21)
        eliJakp(p) = 1;
        eliJak(IDvec) = 1;
    end
end


% Choose observation times to include
inclTime = [0;4;8;12;18;24;36;48;60];

% Storage
relJakMeansIFN = zeros(size(inclTime));
relJakMeansIFNCI = zeros(size(relJakMeansIFN,1),2);
relJakMeansIFNCIstd = zeros(size(relJakMeansIFN,1),2);
nIFN = zeros(size(relJakMeansIFN));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(relJakMeansIFN)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*IFNID);
    relJakMeansIFN(i) = mean(relJak(ID),'omitnan');
    temp2 = relJak(ID);
    n = length(temp2(~isnan(temp2)));
    nIFN(i) = n;
    tempCI = bootci(nboot,{meanfun,relJak(ID)},'Alpha',0.05);
    relJakMeansIFNCI(i,:) = tempCI';
    tempCI = [relJakMeansIFN(i)-tinv(0.975,n-1)*std(relJak(ID),'omitnan')/sqrt(n);...
              relJakMeansIFN(i)+tinv(0.975,n-1)*std(relJak(ID),'omitnan')/sqrt(n)];
    relJakMeansIFNCIstd(i,:) = tempCI';
end


% Plot
figure();
hold on;
plot(inclTime,relJakMeansIFN,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMeansIFNCI(:,1)' fliplr(relJakMeansIFNCI(:,2)')], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
patch([inclTime' fliplr(inclTime')], [relJakMeansIFNCIstd(:,1)' fliplr(relJakMeansIFNCIstd(:,2)')], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMeansIFNCI(:,1),'-k','linewidth',1.5);
plot(inclTime,relJakMeansIFNCI(:,2),'-k','linewidth',1.5);
ylim([0,1.5]);
title('Temporal Evolution of Mean Relative Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
relJakMeansHU = zeros(size(inclTime));
relJakMeansHUCI = zeros(size(relJakMeansHU,1),2);
nHU = zeros(size(relJakMeansHU));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(relJakMeansHU)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*HUID);
    relJakMeansHU(i) = mean(relJak(ID),'omitnan');
    temp2 = relJak(ID);
    n = length(temp2(~isnan(temp2)));
    nHU(i) = n;
%     tempCI = [relJakMeansHU(i)-tinv(0.975,n-1)*std(relJak(ID),'omitnan')/sqrt(n);...
%               relJakMeansHU(i)+tinv(0.975,n-1)*std(relJak(ID),'omitnan')/sqrt(n)];
%     relJakMeansHUCI(i,:) = tempCI';
    tempCI = bootci(nboot,{meanfun,relJak(ID)},'Alpha',0.05);
    relJakMeansHUCI(i,:) = tempCI';
end


% Plot
figure();
hold on;
plot(inclTime,relJakMeansHU,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMeansHUCI(:,1)' fliplr(relJakMeansHUCI(:,2)')], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMeansHUCI(:,1),'-k','linewidth',1.5);
plot(inclTime,relJakMeansHUCI(:,2),'-k','linewidth',1.5);
ylim([0,1.5]);
title('Temporal Evolution of Mean Relative Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage mix
relJakMeansMix = zeros(size(inclTime));
relJakMeansMixCI = zeros(size(relJakMeansMix,1),2);
nMix = zeros(size(relJakMeansMix));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(relJakMeansMix)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*MixID);
    relJakMeansMix(i) = mean(relJak(ID),'omitnan');
    temp2 = relJak(ID);
    n = length(temp2(~isnan(temp2)));
    nMix(i) = n;
%     tempCI = [relJakMeansMix(i)-tinv(0.975,n-1)*std(relJak(ID),'omitnan')/sqrt(n);...
%               relJakMeansMix(i)+tinv(0.975,n-1)*std(relJak(ID),'omitnan')/sqrt(n)];
%     relJakMeansMixCI(i,:) = tempCI';
    tempCI = bootci(nboot,{meanfun,relJak(ID)},'Alpha',0.05);
    relJakMeansMixCI(i,:) = tempCI';
end

% Plot
figure();
hold on;
plot(inclTime,relJakMeansMix,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMeansMixCI(:,1)' fliplr(relJakMeansMixCI(:,2)')], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMeansMixCI(:,1),'-k','linewidth',1.5);
plot(inclTime,relJakMeansMixCI(:,2),'-k','linewidth',1.5);
ylim([0,1.5]);
title('Temporal Evolution of Mean Relative Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('$\frac{JAK2}{JAK2_0}$','fontsize',16,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

% Storage
relJakMeansPV = zeros(size(inclTime));
relJakMeansPVCI = zeros(size(relJakMeansPV,1),2);
relJakMeansPVsd = zeros(size(relJakMeansPV,1),1);
relJakMeansET = zeros(size(inclTime));
relJakMeansETCI = zeros(size(relJakMeansET,1),2);
relJakMeansETsd = zeros(size(relJakMeansET,1),1);
relJakMeansMF = zeros(size(inclTime));
relJakMeansMFCI = zeros(size(relJakMeansMF,1),2);
relJakMeansMFsd = zeros(size(relJakMeansMF,1),1);

% Parameters for bootstrap
nboot = 10000;
medfun = @(x) mean(x,'omitnan');

% Calculate means at different times and confidence intervals
for i=1:length(relJakMeansPV)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp.*PVID.*IFNID);
    relJakMeansPV(i) = mean(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMeansPVsd(i) = std(relJak(ID),'omitnan');
    relJakMeansPVCI(i,:) = tempCI';
    ID = logical(temp.*ETID.*IFNID);
    relJakMeansET(i) = mean(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMeansETCI(i,:) = tempCI';
    relJakMeansETsd(i) = std(relJak(ID),'omitnan');
    ID = logical(temp.*MFID.*IFNID);
    relJakMeansMF(i) = mean(relJak(ID),'omitnan');
    tempCI = bootci(nboot,{medfun,relJak(ID)},'Alpha',0.05);
    relJakMeansMFCI(i,:) = tempCI';
    relJakMeansMFsd(i) = std(relJak(ID),'omitnan');
end

% Plot
figure();
hold on;
plot(inclTime,relJakMeansPV,'.-r','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMeansPVCI(:,1)' fliplr(relJakMeansPVCI(:,2)')], ...
      'r', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMeansPVCI(:,1),'-r','linewidth',1.5);
plot(inclTime,relJakMeansPVCI(:,2),'-r','linewidth',1.5);
plot(inclTime,relJakMeansET,'.-g','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMeansETCI(:,1)' fliplr(relJakMeansETCI(:,2)')], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMeansETCI(:,1),'-g','linewidth',1.5);
plot(inclTime,relJakMeansETCI(:,2),'-g','linewidth',1.5);
plot(inclTime,relJakMeansMF,'.-b','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [relJakMeansMFCI(:,1)' fliplr(relJakMeansMFCI(:,2)')], ...
      'b', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,relJakMeansMFCI(:,1),'-b','linewidth',1.5);
plot(inclTime,relJakMeansMFCI(:,2),'-b','linewidth',1.5);
ylim([0,1.5]);
title('Mean Relative Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('$\frac{JAK2}{JAK2_0}$','fontsize',16,'interpreter','latex');
legend({'PV - Sample relative mean','PV - Bootstrap CI','','',...
        'ET - Sample relative mean','ET - Bootstrap CI','','',...
        'MF - Sample relative mean','MF - Bootstrap CI','',''},...
        'fontsize',12,'location','best','interpreter','latex');

figure();
hold on;
errorbar(inclTime,relJakMeansPV,relJakMeansPVsd,'.-r','linewidth',1.5,'markersize',20);
errorbar(inclTime,relJakMeansET,relJakMeansETsd,'.-g','linewidth',1.5,'markersize',20);
errorbar(inclTime,relJakMeansMF,relJakMeansMFsd,'.-b','linewidth',1.5,'markersize',20);
ylim([0,1.5]);
title('Mean Allele Burden','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean allele burden/\%','fontsize',12,'interpreter','latex');
legend({'PV - Sample relative mean $\pm$ standard deviation'
        'ET - Sample relative mean $\pm$ standard deviation'
        'MF - Sample relative mean $\pm$ standard deviation',},...
        'fontsize',12,'location','best','interpreter','latex');

%% Patient plots - single patient
% Define ID for the patient
pID = 24;
IDvec = patientID==pID;

% Find the baseline date and hospital ID
basDate = dates(logical(IDvec.*baselineID));
hospIDp = hospID(pID);

% Time vectors for plotting
time = relTime(IDvec);

% Vectors for normal regions
timeLin = linspace(minRelTime,maxRelTime+5,2);
minHCTm = 0.40; % fraction, sundhed.dk
maxHCTm = 0.50; % fraction, sundhed.dk
minHCTf = 0.35; % fraction, sundhed.dk
maxHCTf = 0.44; % fraction, sundhed.dk
minHCTmLin = linspace(minHCTm,minHCTm,2);
maxHCTmLin = linspace(maxHCTm,maxHCTm,2);
minHCTfLin = linspace(minHCTf,minHCTf,2);
maxHCTfLin = linspace(maxHCTf,maxHCTf,2);
minHGBm = 8.3; % mmol/L, sundhed.dk
maxHGBm = 10.5; % mmol/L, sundhed.dk
minHGBf = 7.3; % mmol/L, sundhed.dk
maxHGBf = 9.5; % mmol/L, sundhed.dk
minHGBmLin = linspace(minHGBm,minHGBm,2);
maxHGBmLin = linspace(maxHGBm,maxHGBm,2);
minHGBfLin = linspace(minHGBf,minHGBf,2);
maxHGBfLin = linspace(maxHGBf,maxHGBf,2);
minLeuko = 3.5; % *10^9/L, sundhed.dk
maxLeuko = 8.8; % *10^9/L, sundhed.dk
minLeukoLin = linspace(minLeuko,minLeuko,2);
maxLeukoLin = linspace(maxLeuko,maxLeuko,2);
minTRCm = 145; % *10^9/L, sundhed.dk
maxTRCm = 350; % *10^9/L, sundhed.dk
minTRCf = 165; % *10^9/L, sundhed.dk
maxTRCf = 390; % *10^9/L, sundhed.dk
minTRCmLin = linspace(minTRCm,minTRCm,2);
maxTRCmLin = linspace(maxTRCm,maxTRCm,2);
minTRCfLin = linspace(minTRCf,minTRCf,2);
maxTRCfLin = linspace(maxTRCf,maxTRCf,2);
minNeu = 2.0; % *10^9/L, sundhed.dk
maxNeu = 8.8; % *10^9/L, sundhed.dk
minNeuLin = linspace(minNeu,minNeu,2);
maxNeuLin = linspace(maxNeu,maxNeu,2);

% Plot
% Make full screen figure
fh = figure();
fh.WindowState = 'maximized';
tiledlayout(4,3);

nexttile;
stairs(time,IFNsys(IDvec),'-r','linewidth',1.5);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('Dose/$\mu$g average per day','fontsize',14,'interpreter','latex');
title('IFN Dose - Pegasys','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(max(IFNsys),max(IFNintron))]);

nexttile;
stairs(time,IFNintron(IDvec),'-r','linewidth',1.5);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('Dose/$\mu$g average per day','fontsize',14,'interpreter','latex');
title('IFN Dose - PegIntron','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(max(IFNsys),max(IFNintron))]);

nexttile;
stairs(time,HU(IDvec),'-r','linewidth',1.5);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('Dose/mg per day','fontsize',14,'interpreter','latex');
title('HU Dose - Hydrea','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(HU)]);

nexttile;
plot(time,HCT(IDvec),'.r','markersize',20);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('HCT/1','fontsize',14,'interpreter','latex');
title('Haematocrit','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(HCT)]);
if malep(pID)
    patch([timeLin fliplr(timeLin)], [minHCTmLin fliplr(maxHCTmLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
else
    patch([timeLin fliplr(timeLin)], [minHCTfLin fliplr(maxHCTfLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
end

nexttile;
plot(time,HGB(IDvec),'.r','markersize',20);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('HGB/mmol per L','fontsize',14,'interpreter','latex');
title('Haemoglobin','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(HGB)]);
if malep(pID)
    patch([timeLin fliplr(timeLin)], [minHGBmLin fliplr(maxHGBmLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
else
    patch([timeLin fliplr(timeLin)], [minHGBfLin fliplr(maxHGBfLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
end

nexttile;
plot(time,WBC(IDvec),'.r','markersize',20);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('WBC/$10^9$ per L','fontsize',14,'interpreter','latex');
title('Leukocytes','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(WBC)]);
patch([timeLin fliplr(timeLin)], [minLeukoLin fliplr(maxLeukoLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');

nexttile;
plot(time,TRC(IDvec),'.r','markersize',20);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('TRC/$10^9$ per L','fontsize',14,'interpreter','latex');
title('Thrombocytes','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(TRC)]);
if malep(pID)
    patch([timeLin fliplr(timeLin)], [minTRCmLin fliplr(maxTRCmLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
else
    patch([timeLin fliplr(timeLin)], [minTRCfLin fliplr(maxTRCfLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
end

nexttile;
plot(time,neu(IDvec),'.r','markersize',20);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('Neutrophils/$10^9$ per L','fontsize',14,'interpreter','latex');
title('Neutrophils','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(neu)]);
patch([timeLin fliplr(timeLin)], [minNeuLin fliplr(maxNeuLin)], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');

nexttile;
plot(time,jak(IDvec)*100,'.r','markersize',20);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('JAK2/\%','fontsize',14,'interpreter','latex');
title('JAK2 Allele Burden','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,100]);

nexttile;
stairs(time,phlebot(IDvec),'-r','linewidth',1.5);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('Phlebotomies/1','fontsize',14,'interpreter','latex');
title('Phlebotomies','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(phlebot)]);

nexttile;
stairs(time,RBCTrans(IDvec),'-r','linewidth',1.5);
xlabel('Time/months','fontsize',14,'interpreter','latex');
ylabel('RBC Transfusions/1','fontsize',14,'interpreter','latex');
title('Red Blood Cell Transfusions','fontsize',18,'interpreter','latex');
xlim([minRelTime,maxRelTime]);
ylim([0,1.2*max(RBCTrans)]);

%% Patient plots - all patients
% Vectors for normal regions
timeLin = linspace(minRelTime,maxRelTime+5,2);
minHCTm = 0.40; % fraction, sundhed.dk
maxHCTm = 0.50; % fraction, sundhed.dk
minHCTf = 0.35; % fraction, sundhed.dk
maxHCTf = 0.44; % fraction, sundhed.dk
minHCTmLin = linspace(minHCTm,minHCTm,2);
maxHCTmLin = linspace(maxHCTm,maxHCTm,2);
minHCTfLin = linspace(minHCTf,minHCTf,2);
maxHCTfLin = linspace(maxHCTf,maxHCTf,2);
minHGBm = 8.3; % mmol/L, sundhed.dk
maxHGBm = 10.5; % mmol/L, sundhed.dk
minHGBf = 7.3; % mmol/L, sundhed.dk
maxHGBf = 9.5; % mmol/L, sundhed.dk
minHGBmLin = linspace(minHGBm,minHGBm,2);
maxHGBmLin = linspace(maxHGBm,maxHGBm,2);
minHGBfLin = linspace(minHGBf,minHGBf,2);
maxHGBfLin = linspace(maxHGBf,maxHGBf,2);
minLeuko = 3.5; % *10^9/L, sundhed.dk
maxLeuko = 8.8; % *10^9/L, sundhed.dk
minLeukoLin = linspace(minLeuko,minLeuko,2);
maxLeukoLin = linspace(maxLeuko,maxLeuko,2);
minTRCm = 145; % *10^9/L, sundhed.dk
maxTRCm = 350; % *10^9/L, sundhed.dk
minTRCf = 165; % *10^9/L, sundhed.dk
maxTRCf = 390; % *10^9/L, sundhed.dk
minTRCmLin = linspace(minTRCm,minTRCm,2);
maxTRCmLin = linspace(maxTRCm,maxTRCm,2);
minTRCfLin = linspace(minTRCf,minTRCf,2);
maxTRCfLin = linspace(maxTRCf,maxTRCf,2);
minNeu = 2.0; % *10^9/L, sundhed.dk
maxNeu = 8.8; % *10^9/L, sundhed.dk
minNeuLin = linspace(minNeu,minNeu,2);
maxNeuLin = linspace(maxNeu,maxNeu,2);

% Make full screen figure
fh = figure();
fh.WindowState = 'maximized';

% Plot in loop
for p=1:P
    % Define ID for the patient
    pID = p;
    IDvec = patientID==pID;
    
    % Find the baseline date and hospital ID
    basDate = dates(logical(IDvec.*baselineID));
    hospIDp = hospID(pID);
    
    % Time vectors for plotting
    time = relTime(IDvec);

    % Clear figure before plotting
    clf;
    
    % Plot
    tiledlayout(4,3);

    nexttile;
    stairs(time,IFNsys(IDvec),'-r','linewidth',1.5);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('Dose/$\mu$g average per day','fontsize',14,'interpreter','latex');
    title('IFN Dose - Pegasys','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(max(IFNsys),max(IFNintron))]);
    
    nexttile;
    stairs(time,IFNintron(IDvec),'-r','linewidth',1.5);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('Dose/$\mu$g average per day','fontsize',14,'interpreter','latex');
    title('IFN Dose - PegIntron','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(max(IFNsys),max(IFNintron))]);
    
    nexttile;
    stairs(time,HU(IDvec),'-r','linewidth',1.5);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('Dose/mg per day','fontsize',14,'interpreter','latex');
    title('HU Dose - Hydrea','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(HU)]);
    
    nexttile;
    plot(time,HCT(IDvec),'.r','markersize',20);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('HCT/1','fontsize',14,'interpreter','latex');
    title('Haematocrit','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(HCT)]);
    if malep(pID)
        patch([timeLin fliplr(timeLin)], [minHCTmLin fliplr(maxHCTmLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    else
        patch([timeLin fliplr(timeLin)], [minHCTfLin fliplr(maxHCTfLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    end
    
    nexttile;
    plot(time,HGB(IDvec),'.r','markersize',20);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('HGB/mmol per L','fontsize',14,'interpreter','latex');
    title('Haemoglobin','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(HGB)]);
    if malep(pID)
        patch([timeLin fliplr(timeLin)], [minHGBmLin fliplr(maxHGBmLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    else
        patch([timeLin fliplr(timeLin)], [minHGBfLin fliplr(maxHGBfLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    end
    
    nexttile;
    plot(time,WBC(IDvec),'.r','markersize',20);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('WBC/$10^9$ per L','fontsize',14,'interpreter','latex');
    title('Leukocytes','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(WBC)]);
    patch([timeLin fliplr(timeLin)], [minLeukoLin fliplr(maxLeukoLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    
    nexttile;
    plot(time,TRC(IDvec),'.r','markersize',20);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('TRC/$10^9$ per L','fontsize',14,'interpreter','latex');
    title('Thrombocytes','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(TRC)]);
    if malep(pID)
        patch([timeLin fliplr(timeLin)], [minTRCmLin fliplr(maxTRCmLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    else
        patch([timeLin fliplr(timeLin)], [minTRCfLin fliplr(maxTRCfLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    end
    
    nexttile;
    plot(time,neu(IDvec),'.r','markersize',20);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('Neutrophils/$10^9$ per L','fontsize',14,'interpreter','latex');
    title('Neutrophils','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(neu)]);
    patch([timeLin fliplr(timeLin)], [minNeuLin fliplr(maxNeuLin)], ...
          'g', 'FaceAlpha',0.3, 'EdgeColor','none');
    
    nexttile;
    plot(time,jak(IDvec)*100,'.r','markersize',20);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('JAK2/\%','fontsize',14,'interpreter','latex');
    title('JAK2 Allele Burden','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,100]);

    nexttile;
    stairs(time,phlebot(IDvec),'-r','linewidth',1.5);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('Phlebotomies/1','fontsize',14,'interpreter','latex');
    title('Phlebotomies','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(phlebot)]);
    
    nexttile;
    stairs(time,RBCTrans(IDvec),'-r','linewidth',1.5);
    xlabel('Time/months','fontsize',14,'interpreter','latex');
    ylabel('RBC Transfusions/1','fontsize',14,'interpreter','latex');
    title('Red Blood Cell Transfusions','fontsize',18,'interpreter','latex');
    xlim([minRelTime,maxRelTime]);
    ylim([0,1.2*max(RBCTrans)]);
    
    % Save figure
    % saveas(gcf,fullfile('H:\Kode\Billeder og video\DALIAH_5y - Patient Plots\',...
           % sprintf('p%0d_DALIAH_5y.jpg',p)),'jpg');
end

% Close figure
close;

%% IFN means
% Choose observation times to include
inclTime = [0;4;8;12;24];

% Storage
IFNMeans = zeros(size(inclTime));
IFNMeansCI = zeros(size(IFNMeans,1),2);
IFNMeansCIstd = zeros(size(IFNMeans,1),2);
nEli = zeros(size(IFNMeans));

% Parameters for bootstrap
nboot = 10000;
meanfun = @(x) mean(x,'omitnan');

% Calculate means at different times and standard deviations
for i=1:length(IFNMeans)
    temp = timeStudyVisit == inclTime(i);
    ID = logical(temp);
    IFNMeans(i) = mean(IFN(ID),'omitnan');
    temp2 = IFN(ID);
    n = length(temp2(~isnan(temp2)));
    nEli(i) = n;
    tempCI = bootci(nboot,{meanfun,IFN(ID)},'Alpha',0.05);
    IFNMeansCI(i,:) = tempCI';
    tempCI = [IFNMeans(i)-tinv(0.975,n-1)*std(IFN(ID),'omitnan')/sqrt(n);...
              IFNMeans(i)+tinv(0.975,n-1)*std(IFN(ID),'omitnan')/sqrt(n)];
    IFNMeansCIstd(i,:) = tempCI';
end


% Plot
figure();
hold on;
plot(inclTime,IFNMeans,'.-k','linewidth',1.5,'markersize',20);
patch([inclTime' fliplr(inclTime')], [IFNMeansCI(:,1)' fliplr(IFNMeansCI(:,2)')], ...
      'k', 'FaceAlpha',0.3, 'EdgeColor','none');
patch([inclTime' fliplr(inclTime')], [IFNMeansCIstd(:,1)' fliplr(IFNMeansCIstd(:,2)')], ...
      'g', 'FaceAlpha',0.3, 'EdgeColor','none');
plot(inclTime,IFNMeansCI(:,1),'-k','linewidth',1.5);
plot(inclTime,IFNMeansCI(:,2),'-k','linewidth',1.5);
ylim([0,inf]);
title('Temporal Evolution of Mean IFN','fontsize',16,...
      'interpreter','latex');
xlabel('Time/months','fontsize',12,'interpreter','latex');
ylabel('Mean IFN/unit','fontsize',12,'interpreter','latex');
legend({'Sample mean','95\% CI'},'fontsize',12,'location','best',...
        'interpreter','latex');

%% Save data
save('M:\MATLAB\DALIAH_5y_data');