%% Load data
% Change path and save old path
oldFolder = cd('M:\MATLAB');

% Load data
load DALIAH_5y_data;

% Change path back
cd('C:\Users\boklund\Documents\Kode');

%% Add extraTreatIFN, extraTreatRUX, and extraTreatDates
extraTreatDates = cell(P,1);
extraTreatIFN = cell(P,1);
extraTreatRUX = cell(P,1);

%% Export for Julia
% ID for patients we want to export - IFN and at least 3 observations
ID = logical(IFNID.*eliID.*~isnan(relTime));

% Create table
dataTable = table;
dataTable.patientID = patientID(ID);
dataTable.days = relTime(ID)*30.4;
dataTable.IFN = IFN(ID);
dataTable.RUX = zeros(size(IFN(ID)));
dataTable.JAK = jak(ID);
dataTable.TRC = TRC(ID);
dataTable.WBC = WBC(ID);
dataTable.studyVisit = timeStudyVisit(ID);
writetable(dataTable,'M:\data_cancer\DALIAH\Trine_DALIAH_5y\DALIAH_5y - Export for Julia.csv')

%% Numbers for IFN patients with prior HU treatment
% Extract ID for patients fitted
fitIDp = logical(pEliID.*IFNIDp);

% Extract ID for patients fitted and treated with HU prior
priorHUFitIDp = logical(fitIDp.*priorHU);

% Unique patients
unique_patients = unique(patientID);

% Show patients who have received prior HU treatment in the fit
unique_patients(priorHUFitIDp)

%% Check stats of relevant patients
chosenID = false(P,1);
chosenID([34,61,73,80,84,87,187,193,194]) = true;

% age(chosenID)
malep(chosenID)

%% Check number of phlebotomies
% Storage
numPhlebot = zeros(P,1);
for i=1:P
    pID = i;
    tempID = patientID == pID;
    numPhlebot(i) = sum(phlebot(tempID),'omitnan');
end

%% Create table with patient information
% Create table
dataTable2 = table;
dataTable2.patientID = unique_patients(fitIDp);
dataTable2.male = malep(fitIDp);
dataTable2.age = age(fitIDp);
dataTable2.PV = PVIDp(fitIDp);
dataTable2.ET = ETIDp(fitIDp);
dataTable2.MF = MFIDp(fitIDp);
dataTable2.Phlebotomies = numPhlebot(fitIDp);
dataTable2.initVAF = basJak(fitIDp);
dataTable2.Pegasys = sysIDp(fitIDp);
dataTable2.PegIntron = intronIDp(fitIDp);
dataTable2.BothIFN = sysIntronIDp(fitIDp);

% Save table
writetable(dataTable2,'M:\data_cancer\DALIAH\Trine_DALIAH_5y\DALIAH_5y_patient_information - Export for Julia.csv')