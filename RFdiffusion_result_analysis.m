function RFdiffusion_result_analysis
disp('           ')
disp('Program start. BaeLab RFdiffuion analysis matlab script. Contact Heesoo for any issue.')
%% File directories and general settings
OriginalDirectory=cd;

directory_path{1} = 'C:\PE-SB_example\PE-SB_ExperimentCandidates\MCA';

% directory_path{1} = 'C:\PE-SB_example\PE-SB_ExperimentCandidates\MC';

% directory_path{1} = 'C:\PE-SB_example\PE-SB_TotalCandidates\MCA';

% directory_path{1} = 'C:\PE-SB_example\PE-SB_TotalCandidates\MC';

rmsd_limit = 1.5;
reverse_rmsd_limit = 'n';
plddt_limit = 0.85;
i_pae_limit = 8.5;

peptidelength_limit = [50, 150];

topranknumber =23;
toprank_save = 'n';             %to save top rank files seperately
toprank_savepath = 'C:\PE-SB_example';
toprank_savehead = 'MCA';

sortby = 'i';                   %rank sort by 'r' for rmsd or 'i' for i_pae

usemolviewer = 'n';             %to see structures
peptide_or_competition = 'p';   %for molviewer, 'p' : peptide binding dimer structure, 'c' : competition trimer structure, 'd' : designed peptide binding dimer structure
savejsonfile = 'n';             %to generate AF3 input files
useAF3competitionfiles = 'y';   %when you have AF3 result files

ranktext_on = 'y';              %to see rank text on figures
RMSD_onAF3 = 'y';               %to see RMSD color bar on AF3 figure
filepath_on = 'n';              %to see filepath of top ranks


%% AF3 competition settings
if useAF3competitionfiles == 'y'
    AF3_bindercompetition_path{1} = 'C:\PE-SB_example\PE-SB_ExperimentCandidates\MCA\MCA_AF3data';
    % AF3_bindercompetition_path{1} = 'C:\PE-SB_example\PE-SB_ExperimentCandidates\MC\MC_AF3data';
    % AF3_bindercompetition_path{1} = 'C:\PE-SB_example\PE-SB_TotalCandidates\MCA\MCA_AF3data';
    % AF3_bindercompetition_path{1} = 'C:\PE-SB_example\PE-SB_TotalCandidates\MC\MC_AF3data';
    
    AtomPairs = [2191 2995; 2047 4230; 436 2857; 2080 4176; 449 2818]   %to specify atom pairs. Careful, these are not residue numbers, it's atom numbers in cif file(=AF3 output).
    
    pairnumber = length(AtomPairs(:,1));

    AF3originalcomplex_filepath = 'C:\PE-SB_example\PE-SB_ExperimentCandidates\MLH1CPMS2C.cif'    %original dimer structure
    AtomPairDistanceOriginal = AtomPairDistance(AF3originalcomplex_filepath, AtomPairs, pairnumber)
    OriginalDistanceAverage = mean(AtomPairDistanceOriginal)

    distancelimit = 10^0.99;      %to define successful inhibition
    reverse_distancelimit = 'n';   %to investigate failed candidates, use 'y'
end



%% AF3 json file settings
if savejsonfile == 'y'
    output_path = 'C:\PE-SB_example';

    jsonnamehead = 'MCA_MLH1C_binder';    %to specify job and file name

    targetproteinsequence = 'SRKEMTAACTPRRRIINLTSVLSLQEEINEQGHEVLREMLHNHSFVGCVNPQWALAQHQTKLYLLNTTKLSEELFYQILIYDFANFGVLRLSEPAPLFDLAMLALDSPESGWTEEDGPKEGLAEYIVEFLKKKAEMLADYFSLEIDEEGNLIGLPLLIDNYVPPLEGLPIFILRLATEVNWDEEKECFESLSKECAMFYSIRKQYISEESTLSGQQSEVPGSIPNSWKWTVEHIVYKALRSHILPPKHFTEDGNILQLANLPDLYKVFERC';
    partnerproteinsequence = 'DVAVKINKKVVPLDFSMSSLAKRIKQLHHEAQQSEGEQNYRKFRAKICPGENQAAEDELRKEISKTMFAEMEIIGQFNLGFIITKLNEDIFIVDQHATDEKYNFEMLQQHTVLQGQRLIAPQTLNLTAVNEAVLIENLEIFRKNGFDFVIDENAPVTERAKLISLPTSKNWTFGPQDVDELIFMLSDSPGVMCRPSRVKQMFASRACRKSVMIGTALNTSEMKKLITHMGEMDHPWNCPHGRPTMRHIANLGVISQN';
    binderproteinsequence = 'SSMEERLKRIAERARKIREKEEEEEKKEALKLAENFKADLEKMKKAYKEAVKKGKEVLKKEGYEAALEEVKKVVDEMLKGNKIGQMVVLMLVERELEEIKAREE';     % will be overwritten. this is just an example text to be seen.

    jsoncontent.name = jsonnamehead;
    jsoncontent.modelSeeds = {'2096432081'};    % you can change the seed as prefered
    jsoncontent.sequences{1}.proteinChain.sequence = targetproteinsequence;
    jsoncontent.sequences{1}.proteinChain.count = int16(1);
    jsoncontent.sequences{2}.proteinChain.sequence = partnerproteinsequence;
    jsoncontent.sequences{2}.proteinChain.count = int16(1);
    jsoncontent.sequences{3}.proteinChain.sequence = binderproteinsequence;
    jsoncontent.sequences{3}.proteinChain.count = int16(1);
end
%% Load binder files

disp(['rmsd_limit = ' num2str(rmsd_limit) '   plddt_limit = ' num2str(plddt_limit) '   i_pae_limit = ' num2str(i_pae_limit) '   topranknumber = ' num2str(topranknumber) ]);

NumberofFiles = 0;

directory_number = length(directory_path);

for i = 1:directory_number
    cd(directory_path{i});
    MyFolderInfo = dir('**/*best.pdb');
    filenumber = length(MyFolderInfo);
    for j = 1:filenumber
        NumberofFiles = NumberofFiles + 1;
        pdbfilepath{NumberofFiles} = fullfile(MyFolderInfo(j).folder, MyFolderInfo(j).name);
    end
end

cd(OriginalDirectory);

for i = 1:NumberofFiles
    fid = fopen(pdbfilepath{i});
    line_ex = fgetl(fid);
    fclose(fid);
    texts = split(line_ex);
    file_rmsd_data(i) = str2num(texts{8});
end

NumberofDesign = 0;
for j = 1:NumberofFiles
    FileInfo = dir(pdbfilepath{j});
    designfilepath = [pdbfilepath{j}(1:(end-8)) 'design.fasta'];
    [fastaheader, fastasequence] = fastaread(designfilepath);
    if iscell(fastaheader) == 0
        temp = fastaheader;
        clear fastaheader;
        fastaheader{1} = temp;
        temp = fastasequence;
        clear fastasequence;
        fastasequence{1} = temp;
    end
    for i = 1:length(fastaheader)
        headerarray = split(fastaheader{i},[":", "|"]);
        if file_rmsd_data(j) == str2num(headerarray{end})
            sequences= split(fastasequence{i},'/');
            if length(sequences{2}) >= peptidelength_limit(1) && length(sequences{2}) <= peptidelength_limit(2)
                NumberofDesign = NumberofDesign + 1;
                pdb_filepath{NumberofDesign} = pdbfilepath{j};
                design_filepath{NumberofDesign} = designfilepath;
                rmsd_data(NumberofDesign) = file_rmsd_data(j);
                targetsequences{NumberofDesign} = sequences{1};
                peptidesequences{NumberofDesign} = sequences{2};
                targetlengths(NumberofDesign) = length(sequences{1});
                peptidelengths(NumberofDesign) = length(sequences{2});
                plddts(NumberofDesign) = str2num(headerarray{7});
                i_paes(NumberofDesign) = str2num(headerarray{11});
                break;
            end
        end
    end
end
NumberofFiles
NumberofDesign
%% Load AF3 files
if useAF3competitionfiles == 'y'
    peptide_atompair_distance = zeros(NumberofDesign, pairnumber);
    peptide_atompair_distance_average = zeros(1, NumberofDesign);
    Number_AF3_Files = 0;
    
    AF3_directory_number = length(AF3_bindercompetition_path);
    for i = 1:AF3_directory_number
        cd(AF3_bindercompetition_path{i});
        MyFolderInfo = dir('**/*model_0.cif');
        filenumber = length(MyFolderInfo);
        for j = 1:filenumber
            Number_AF3_Files = Number_AF3_Files + 1;
            AF3filepath{Number_AF3_Files} = fullfile(MyFolderInfo(j).folder, MyFolderInfo(j).name);
        end
    end
    Number_AF3_Files
    
    cd(OriginalDirectory);
    
    for i = 1:Number_AF3_Files
        AtomPair_distance(i,:) = AtomPairDistance(AF3filepath{i}, AtomPairs, pairnumber);
        AtomPair_distanceAverage(i) = mean(AtomPair_distance(i,:));
        jsonfilepath = [AF3filepath{i}(1:(end-11)) 'job_request.json'];
        AF3json = readstruct(jsonfilepath);
        AF3_bindersequence(i) = AF3json.sequences(3).proteinChain.sequence;
        
        binderindex = find(contains(peptidesequences,AF3_bindersequence(i)));
        if ~isempty(binderindex)
            peptide_atompair_distance(binderindex,:) = AtomPair_distance(i,:);
            peptide_atompair_distance_average(binderindex) = AtomPair_distanceAverage(i);
            peptide_AF3filepath{binderindex} = AF3filepath{i};
            peptide_AF3jsonfilepath{binderindex} = jsonfilepath;
        else
            disp(['AF3only without candidate (check peptidelength : ' num2str(length(AF3_bindersequence{i})) ') : ' AF3_bindersequence{i}])
        end
    
    end
    
    AF3_index = (peptide_atompair_distance_average ~= 0);
end

%% Histograms, distributions, and top ranks

hf = figure(100);
hf.Name = 'Histograms for Total candidates';
hf.Color = "white";
subplot(2,2,1);
edges = [0:1:200];
h = histogram(peptidelengths,edges);
title('Peptide length histogram')
xlabel('Length (aa)')
ylabel('Count')
axis([0 200 0 inf])

subplot(2,2,2);
edges = [0:2:100];
h = histogram(rmsd_data,edges);
title('RMSD histogram')
xlabel('RMSD (Å)')
ylabel('Count')
axis([0 100 0 inf])

subplot(2,2,3);
edges = [0:0.02:1];
h = histogram(plddts,edges);
title('pLDDT histogram')
xlabel('pLDDT')
ylabel('Count')
axis([0 1 0 inf])

subplot(2,2,4);
edges = [0:1:50];
h = histogram(i_paes,edges);
title('interaction pAE histogram')
xlabel('interaction pAE')
ylabel('Count')
axis([0 50 0 inf])



total_number = length(plddts);

hf = figure(101);
hf.Name = 'pLDDT i_pAE RMSD filtering';
hf.Color = "white";
hf.Position(3:4) = [1200 400];
subplot(1,2,1);
sz = 25;
scatter(plddts, i_paes, sz, rmsd_data, "filled");
title(['Total number = ' num2str(total_number)])
xlabel('pLDDT')
ylabel('interaction pAE')
c = colorbar('FontSize', 14);
c.Label.String = 'RMSD (Å)';
axis([0.3 1 0 35])
ax = gca; ax.LineWidth = 2; ax.FontSize  = 14;
clim([0 5])
boxright = 0.95;
boxbottom = 5;
line([plddt_limit plddt_limit boxright boxright plddt_limit],[boxbottom i_pae_limit i_pae_limit boxbottom boxbottom], 'LineWidth', 2)

if reverse_rmsd_limit ~= 'y'
    good_candidate_index = logical((plddts>=plddt_limit ).*(i_paes<=i_pae_limit).*(rmsd_data<=rmsd_limit));
else
    good_candidate_index = logical((plddts>=plddt_limit ).*(i_paes<=i_pae_limit).*(rmsd_data>rmsd_limit));
end
good_candidate_number = sum(good_candidate_index);
good_candidate_percent = good_candidate_number/total_number;


subplot(1,2,2);
sz = 50;
scatter(plddts(good_candidate_index), i_paes(good_candidate_index), sz, rmsd_data(good_candidate_index), "filled");
title(['Filtered Candidates = ' num2str(good_candidate_number) '/' num2str(total_number) '  ratio = ' num2str(good_candidate_percent)])
xlabel('pLDDT')
ylabel('interaction pAE')
c = colorbar('FontSize', 14);
c.Label.String = 'RMSD (Å)';
axis([plddt_limit boxright boxbottom i_pae_limit])
ax = gca; ax.LineWidth = 2; ax.FontSize  = 14;
clim([0 rmsd_limit])
% grid on
% grid minor



if useAF3competitionfiles == 'y'
    MissingStructure_index = AF3_index - good_candidate_index;
    NoAF3_in_goodcandidates = find(MissingStructure_index == -1);
    AF3_in_badcandidates = find(MissingStructure_index == 1);
    number_NoAF3_in_goodcandidates = length(NoAF3_in_goodcandidates);
    for i = 1:number_NoAF3_in_goodcandidates
        if filepath_on == 'y'
            disp(['No AF3 structure (path) : ' pdb_filepath{NoAF3_in_goodcandidates(i)}]);
        end
        disp(['No AF3 structure : ' num2str(i) ' / ' num2str(number_NoAF3_in_goodcandidates) '   Length : ' num2str(targetlengths(NoAF3_in_goodcandidates(i))) ' ' num2str(peptidelengths(NoAF3_in_goodcandidates(i))) '   rmsd: ' num2str(rmsd_data(NoAF3_in_goodcandidates(i))) '  plddt: ' num2str(plddts(NoAF3_in_goodcandidates(i))) '  i_pae: ' num2str(i_paes(NoAF3_in_goodcandidates(i))) '  dist.: ' num2str(peptide_atompair_distance_average(NoAF3_in_goodcandidates(i))) '  seq.: ' peptidesequences{NoAF3_in_goodcandidates(i)} ])
                    
        if savejsonfile == 'y'
            jsoncontent.name = [jsonnamehead ' ' num2str(i) '_noAF3_' num2str(number_NoAF3_in_goodcandidates) ' ' peptidesequences{NoAF3_in_goodcandidates(i)}(1:5)];
            jsoncontent.sequences{3}.proteinChain.sequence = peptidesequences{NoAF3_in_goodcandidates(i)};
            writelines(["[ " jsonencode(jsoncontent, PrettyPrint=true) " ]"], fullfile(output_path, [lower(replace(jsoncontent.name, ["-"; ":"; " "],"_")) '_job_request.json']));
        end
    end
    number_AF3_in_badcandidates = length(AF3_in_badcandidates);
    for i = 1:number_AF3_in_badcandidates
        if filepath_on == 'y'
            disp(['AF3 structure in bad candidates (path) : ' pdb_filepath{AF3_in_badcandidates(i)}]);
        end
        disp(['AF3 structure in bad candidates : ' num2str(i) ' / ' num2str(number_AF3_in_badcandidates) '   Length : ' num2str(targetlengths(AF3_in_badcandidates(i))) ' ' num2str(peptidelengths(AF3_in_badcandidates(i))) '   rmsd: ' num2str(rmsd_data(AF3_in_badcandidates(i))) '  plddt: ' num2str(plddts(AF3_in_badcandidates(i))) '  i_pae: ' num2str(i_paes(AF3_in_badcandidates(i))) '  dist.: ' num2str(peptide_atompair_distance_average(AF3_in_badcandidates(i))) '  seq.: ' peptidesequences{AF3_in_badcandidates(i)} ])
    end

    if reverse_distancelimit ~= 'y'
        AF3_good_candidate_index = (peptide_atompair_distance_average>=distancelimit ).*AF3_index.*good_candidate_index;
    else
        AF3_good_candidate_index = (peptide_atompair_distance_average<distancelimit ).*AF3_index.*good_candidate_index;
    end
    AF3_good_candidate_number = sum(AF3_good_candidate_index);
    good_candidate_index_withAF3 = logical(AF3_index.*good_candidate_index);

    AF3number = sum(sum(AF3_index));

    hf102 = figure(102);
    hf102.Name = 'AF3 competition filtering';
    hf102.Color = "white";
    hf102.Position(3:4) = [1200 600];
    subplot(1,2,1);
    sz = 50;
    scatter(plddts(AF3_index), i_paes(AF3_index), sz, peptide_atompair_distance_average(AF3_index), "filled");
    title(['AF3 Total number = ' num2str(AF3number)   ', AF3 Filtered Candidates = ' num2str(AF3_good_candidate_number) '/' num2str(good_candidate_number) ])
    xlabel('pLDDT')
    ylabel('interaction pAE')
    c = colorbar('FontSize', 14);
    c.Label.String = 'Distance (Å)';
    axis([plddt_limit 1 0 i_pae_limit])
    ax = gca; ax.LineWidth = 2; ax.FontSize  = 14;
    clim([7 10])
    
    
    subplot(4,2,2);
    AF3_log10_distance = log10(peptide_atompair_distance_average(good_candidate_index_withAF3));
    AF3_i_paes = i_paes(good_candidate_index_withAF3);
    leftvalue = log10(7);
    rightvalue = log10(10^1.8);
    edges = [leftvalue:0.01:rightvalue];
    h = histogram(AF3_log10_distance,edges);
    ylabel('Count')
    axis([leftvalue rightvalue 0 inf])
    axAF3plot = gca; axAF3plot.LineWidth = 2; axAF3plot.FontSize  = 14;
    xticks([]);
    line([log10(OriginalDistanceAverage) log10(OriginalDistanceAverage)],[0 max(h.Values)], 'LineWidth', 2)
    line([log10(distancelimit) log10(distancelimit)],[0 max(h.Values)], 'LineWidth', 2)

    h_subplot = subplot(4,2,[4, 6]);
    sz = 50;

    if RMSD_onAF3 == 'y'
        scatter(log10(peptide_atompair_distance_average(good_candidate_index_withAF3)), i_paes(good_candidate_index_withAF3), sz, rmsd_data(good_candidate_index_withAF3), "filled");
        c = colorbar('FontSize', 14);
        c.Label.String = 'RMSD (Å)';
        clim([0 rmsd_limit])
    else
        scatter(log10(peptide_atompair_distance_average(good_candidate_index_withAF3)), i_paes(good_candidate_index_withAF3), sz, "filled");
    end
    xlabel('log_1_0 distance (Å)')
    ylabel('interaction pAE')
    axis([leftvalue rightvalue 5.5 i_pae_limit])
    ax102 = gca; ax102.LineWidth = 2; ax102.FontSize  = 14;
    
    line([log10(OriginalDistanceAverage) log10(OriginalDistanceAverage)],[0 i_pae_limit], 'LineWidth', 2)
    line([log10(distancelimit) log10(distancelimit)],[0 i_pae_limit], 'LineWidth', 2)

    ax = gca;
    axAF3plot.Position(3) = ax.Position(3);
    
    hf103 = figure(103);
    hf103.Name = 'Atom-pair distance distribution';
    hf103.Color = "white";
    hf103.Position(3:4) = [1200 400];
    AtomPairDistanceOriginal_array = ones(AF3number,pairnumber).*AtomPairDistanceOriginal;
    for i = 1:pairnumber
        ax103(1,i) = subplot(2,pairnumber,i);
        sz = 50;
        scatter(log10(peptide_atompair_distance_average(AF3_index)), log10(peptide_atompair_distance(AF3_index,i)), sz, i_paes(AF3_index), "filled");
        xlabel('log_1_0 distance (Å) average')
        ylabel(['log_1_0 distance (Å) ' num2str(i) 'th pair'])
        c = colorbar;
        c.Label.String = 'interaction pAE';
        axis([log10(5) log10(200) log10(5) log10(200)])
        ax = gca; ax.LineWidth = 2; ax.FontSize  = 14;
        clim([5 i_pae_limit])

        line([0 100],[0 100], 'LineWidth', 1);

        ax103(2,i) = subplot(2,pairnumber,pairnumber+i);
        sz = 50;
        scatter(log10(AtomPairDistanceOriginal_array(:,i)), log10(peptide_atompair_distance(AF3_index,i)), sz, i_paes(AF3_index), "filled");
        xlabel('log_1_0 distance (Å) original')
        ylabel(['log_1_0 distance (Å) ' num2str(i) 'th pair'])
        c = colorbar;
        c.Label.String = 'interaction pAE';
        axis([log10(5) log10(200) log10(5) log10(200)])
        ax = gca; ax.LineWidth = 2; ax.FontSize  = 14;
        clim([5 i_pae_limit])

        line([0 100],[0 100], 'LineWidth', 1);
    end
    
end


hf = figure(101);
hf.Color = "white";
subplot(1,2,2);


if sortby == 'r'
    [sorted_rmsd, sorted_index]= sort(rmsd_data);
elseif sortby == 'i'
    [sorted_i_paes, sorted_index]= sort(i_paes);
else
    disp("Error : your sortby is not 'r' or 'i'.");
end
filtered_index = [];
ranknumber = 0;
j = 0;
while ranknumber < topranknumber && j<NumberofDesign
    j = j + 1;
    if good_candidate_index(sorted_index(j)) == 1
        binderok = 'n';
        if useAF3competitionfiles == 'y'
            if AF3_good_candidate_index(sorted_index(j)) == 1
                binderok = 'y';
            end
        else 
            binderok = 'y';
        end
        if binderok == 'y'
            ranknumber = ranknumber + 1;
            
            if ranknumber == 1 || ~strcmp(targetsequences{sorted_index(1)}, targetsequences{sorted_index(j)})
                disp(['tartget: ' targetsequences{sorted_index(j)}]);
            end

            if useAF3competitionfiles == 'y'
                disp(['Rank # : ' num2str(ranknumber) ' / ' num2str(NumberofDesign) '   Length : ' num2str(targetlengths(sorted_index(j))) ' ' num2str(peptidelengths(sorted_index(j))) '   rmsd: ' num2str(rmsd_data(sorted_index(j))) '  plddt: ' num2str(plddts(sorted_index(j))) '  i_pae: ' num2str(i_paes(sorted_index(j))) '  dist.: ' num2str(peptide_atompair_distance_average(sorted_index(j))) '  seq.: ' peptidesequences{sorted_index(j)} ])
                if filepath_on == 'y'
                    disp(peptide_AF3filepath{sorted_index(j)})
                end
                if toprank_save == 'y'
                    copyfile(peptide_AF3filepath{sorted_index(j)}, toprank_savepath);
                    [~,name,ext] = fileparts(peptide_AF3filepath{sorted_index(j)});
                    origianlfilename = [name ext];
                    ciffilename = [toprank_savehead '_rank_' num2str(ranknumber, '%04i') '_trimer_model_0.cif'];
                    if strcmp(origianlfilename, ciffilename) == 0
                        movefile(fullfile(toprank_savepath, origianlfilename), fullfile(toprank_savepath, ciffilename));
                    end
                    AF3json = readstruct(peptide_AF3jsonfilepath{sorted_index(j)});
                    AF3json.name = ciffilename(1:(end-12));
                    if iscell(AF3json.modelSeeds) == 0
                        temp = AF3json.modelSeeds;
                        AF3json.modelSeeds = {temp};
                    end
                    writelines(["[ " jsonencode(AF3json, PrettyPrint=true) " ]"], fullfile(toprank_savepath, [ciffilename(1:(end-11)) 'job_request.json']));
                end
            else
                disp(['Rank # : ' num2str(ranknumber) ' / ' num2str(NumberofDesign) '   Length : ' num2str(targetlengths(sorted_index(j))) ' ' num2str(peptidelengths(sorted_index(j))) '   rmsd: ' num2str(rmsd_data(sorted_index(j))) '  plddt: ' num2str(plddts(sorted_index(j))) '  i_pae: ' num2str(i_paes(sorted_index(j))) '  peptide: ' peptidesequences{sorted_index(j)} ])
            end

            if filepath_on == 'y'
                disp(pdb_filepath{sorted_index(j)})
            end
            if toprank_save == 'y'
                copyfile(pdb_filepath{sorted_index(j)}, toprank_savepath);
                [~,name,ext] = fileparts(pdb_filepath{sorted_index(j)});
                origianlfilename = [name ext];
                pdbfilename = [toprank_savehead '_rank_' num2str(ranknumber, '%04i') '_SB_best.pdb'];
                if strcmp(origianlfilename, pdbfilename) == 0
                    movefile(fullfile(toprank_savepath, origianlfilename), fullfile(toprank_savepath, pdbfilename));
                end
                copyfile(design_filepath{sorted_index(j)}, toprank_savepath);
                [~,name,ext] = fileparts(design_filepath{sorted_index(j)});
                origianlfilename = [name ext];
                fastafilename = [pdbfilename(1:(end-8)) 'design.fasta'];
                if strcmp(origianlfilename, fastafilename) == 0
                    movefile(fullfile(toprank_savepath, origianlfilename), fullfile(toprank_savepath, fastafilename));
                end
            end

            filtered_index(ranknumber) = sorted_index(j);
            if ranktext_on == 'y'
                figure(101);
                subplot(1,2,2);
                text(plddts(sorted_index(j)),i_paes(sorted_index(j)),num2str(ranknumber),'FontSize',14);
            end

            if useAF3competitionfiles == 'y' && ranktext_on == 'y'
                if ranktext_on == 'y'
                    figure(102);
                    subplot(ax102);
                    text(log10(peptide_atompair_distance_average(sorted_index(j))), i_paes(sorted_index(j)), num2str(ranknumber),'FontSize',14);
                    for i = 1:pairnumber
                        figure(103);
                        subplot(ax103(1,i));
                        text(log10(peptide_atompair_distance_average(sorted_index(j))), log10(peptide_atompair_distance(sorted_index(j),i)), num2str(ranknumber),'FontSize',14);
                        subplot(ax103(2,i));
                        text(log10(AtomPairDistanceOriginal_array(1,i)), log10(peptide_atompair_distance(sorted_index(j),i)), num2str(ranknumber),'FontSize',14);
                    end
                end
                axAF3plot.Position(3) = ax102.Position(3);
            end

            if savejsonfile == 'y' && useAF3competitionfiles ~= 'y'
                jsoncontent.name = [jsonnamehead ' ' num2str(ranknumber) '_' num2str(NumberofDesign) ' ' peptidesequences{sorted_index(j)}(1:5)];
                jsoncontent.sequences{3}.proteinChain.sequence = peptidesequences{sorted_index(j)};
                writelines(["[ " jsonencode(jsoncontent, PrettyPrint=true) " ]"], fullfile(output_path, [lower(replace(jsoncontent.name, ["-"; ":"; " "],"_")) '_job_request.json']));
            end
        end
    end
end
if usemolviewer == 'y'
    for i = 1:ranknumber
        if peptide_or_competition == 'p'
            molviewer(pdb_filepath{filtered_index(i)});
        elseif peptide_or_competition == 'c'
            molviewer(peptide_AF3filepath{filtered_index(i)});  %Careful, matlab molviewer has some errors for cif files.
        elseif peptide_or_competition == 'd' % for designed model
            [filepath,name,ext] = fileparts(pdb_filepath{filtered_index(i)})
            cd(filepath);
            cd("../");
            designpdb = dir('*.pdb');
            designpdbpath = fullfile(designpdb.folder, designpdb.name);
            molviewer(designpdbpath);
        else
            disp("Error : your peptide_or_competition is not 'p', 'c' or 'd'.")
        end
    end
end

hf = figure(105);
hf.Name = 'Histograms for filtered candidates';
hf.Color = "white";

subplot(2,2,1);
edges = [0:1:200];
h = histogram(peptidelengths(filtered_index),edges);
title('Peptide length histogram')
xlabel('Length (aa)')
ylabel('Count')
axis([0 200 0 inf])

subplot(2,2,2);
edges = [0:2:100];
h = histogram(rmsd_data(filtered_index),edges);
title('RMSD histogram')
xlabel('RMSD (Å)')
ylabel('Count')
axis([0 100 0 inf])

subplot(2,2,3);
edges = [0:0.02:1];
h = histogram(plddts(filtered_index),edges);
title('pLDDT histogram')
xlabel('pLDDT')
ylabel('Count')
axis([0 1 0 inf])

subplot(2,2,4);
edges = [0:1:50];
h = histogram(i_paes(filtered_index),edges);
title('interaction pAE histogram')
xlabel('interaction pAE')
ylabel('Count')
axis([0 50 0 inf])


disp('Program end.')

cd(OriginalDirectory);

end

function AtomPair_distance = AtomPairDistance(AF3originalcomplex_filepath, AtomPairs, pairnumber)
    textlines = readlines(AF3originalcomplex_filepath,"EmptyLineRule","skip","WhitespaceRule","trimleading");
    
    for j = 1:pairnumber
        atomindex = contains(textlines,['ATOM ' num2str(AtomPairs(j,1)) ' ']);
        atomtexts = textlines(atomindex);
        atomnumber = length(atomtexts);
        for k = 1:atomnumber
            atomcontent = split(atomtexts(k));
            atomxyz(j,1,1:3) = [str2num(atomcontent(11)) str2num(atomcontent(12)) str2num(atomcontent(13))];
        end
    
        
        atomindex = contains(textlines,['ATOM ' num2str(AtomPairs(j,2)) ' ']);
        atomtexts = textlines(atomindex);
        atomnumber = length(atomtexts);
        for k = 1:atomnumber
            atomcontent = split(atomtexts);
            atomxyz(j,2,1:3) = [str2num(atomcontent(11)) str2num(atomcontent(12)) str2num(atomcontent(13))];
        end
    
        for k = 1:3
            diff(k) = atomxyz(j,1,k)-atomxyz(j,2,k);
        end
        AtomPair_distance(j) = norm(diff);
    end
end




