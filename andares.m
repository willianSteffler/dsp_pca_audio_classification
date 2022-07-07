clc;clear all; close all;

reload = 0; %se tem que carregar novamente os dados em cada iteração

% ADp = repelem([200],1)';
% Ap = repelem([0.7],length(ADp))';  %define % de PC
% ATst_p = repelem([0.2],length(ADp))'; %quantidade teste
% Aset=repelem(["augmented"],length(ADp))';
% Aconfusion = repelem([0],length(ADp))';
% Aguess_other = repelem([1],length(ADp))';

ATst_p = repelem([0],1)'; %quantidade teste
ADp = repelem([200],length(ATst_p))';
Ap = repelem([1],length(ATst_p))';  %define % de PC
Aset=repelem(["initial"],length(ATst_p))';
Aconfusion = repelem([0],length(ATst_p))';
Aguess_other = repelem([1],length(ATst_p))';

test_cases = table(ATst_p,Ap,ADp,Aset,Aconfusion,Aguess_other);

for tc_row = 1:size(test_cases,1)
    close all;
    curr_case = test_cases(tc_row,:);

    Tst_p = curr_case.ATst_p; %quantidade teste
    p = curr_case.Ap;  %define % de PC
    Dp = curr_case.ADp;
    set= curr_case.Aset;
    confusion = curr_case.Aconfusion;
    guess_other = curr_case.Aguess_other;

    folder = "experiments/"+ set;
    files = dir(folder);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags); % A structure with extra info.
    % Get only the folder names into a cell array.
    subFolderNames = {subFolders(3:end).name};
    exp_no = length(subFolderNames);
    exp_folder = folder+ "/" + string(exp_no);
    
    if ~exist(exp_folder, 'dir')
        mkdir(exp_folder)
    end

    if reload == 1 || tc_row == 1
        
        %amostras primeiro andar
        PA_files = get_files("audio/"+ set +"/primeiro_andar");
        All_files = PA_files;
        %amostras de teste
        PAt_ex = round(Tst_p * length(PA_files));
        tidx = randperm(length(PA_files),PAt_ex);
        PAt_files = PA_files(tidx);
        
        %amostras de treino
        PA_files = PA_files(setdiff(1:length(PA_files), tidx));
        PA_ex = length(PA_files); %numero de exemplos primeiro andar
        
        %amostras segundo andar
        SA_files = get_files("audio/"+ set +"/segundo_andar");
        All_files = [All_files;SA_files];
        %amostras de teste
        SAt_ex = round(Tst_p * length(SA_files));
        tidx = randperm(length(SA_files),SAt_ex);
        SAt_files = SA_files(tidx);
        SAt_off = PAt_ex + 1;
        
        %amostras de treino
        SA_files = SA_files(setdiff(1:length(SA_files), tidx));
        SA_ex = length(SA_files); %indice do inicio de exemplos segundo andar
        SA_off = PA_ex + 1;
        
        %encontra menor amostra
        ma = menor_amostra(All_files);
    end
    
    %------------------------------------------------TREINO------------------------------------------------
    
    Filepath = [PA_files; SA_files];
    Tex = length(Filepath)
    
    fprintf('###########################################\n')
    X = zeros(ma,Tex);
    %reduz dimensão das amostras
    for i = 1 : Tex
        Filepath{i,1};
        S = resize(Filepath{i,1},ma);
        X(:,i) = S;
    end
    
    %size(X)
    
    mu = mean(X,2);  %Média para cada ponto P de todas as amostras (2)
    
    %wt= V'*(K-mu) Padroes Centralizados (3)
    Phi = X;
    for i = 1:Tex
       Phi(:,i) = Phi(:,i) - mu;
    end
    
    Cxx = (Phi'*Phi);   %matriz de covariancia (8)
    
    %size(Cxx)
    
    [U, Lambda] = eig(Cxx);
    
    %organiza autovalores
    [l, ind] = sort(diag(Lambda),"descend");
    Lambda = Lambda(ind, ind);
    U = U(:, ind);
    U = U(:, 1:Tex-1);
    tr = sum(l);
    
    figure
    bar(l./tr*100)
    title('Principais componentes')
    ylabel('% de variancia') 
    xlabel('PCn') 
    saveas(gcf,exp_folder+'/PC.png')
    
    for f = 1:Tex
      if (sum(l(1:f)) / tr >= p)
        break;
      end
    end
    
    % Remove os autovalores de menor energia
    Lambda = Lambda(:, 1:f);
    Lambda = Lambda(1:f, :);
    
    % Remove os autovetores de menor energia
    U = U(:, 1:f);
    
    %disp(norm(U(:, 1)));
    V = Phi * U; % (10)
    %disp(norm(V(:, 1)));
    %normalização
    for i = 1:f
      V(:, i) = V(:, i) / norm(V(:, i));
    end
    
    %projeção das variaveis nos componmentes principais
    PC3_var = sum(l(1:3))/tr
    PC2_var = sum(l(1:2))/tr
    
    ux_p = project(U,X,1,PA_ex);
    ux_s = project(U,X,SA_off,Tex);
    
    figure('Renderer', 'painters', 'Position', [10 10 650 460])
    scatter3(ux_p(:,1),ux_p(:,2),ux_p(:,3),5,'green')
    hold on;
    scatter3(ux_s(:,1),ux_s(:,2),ux_s(:,3),5,'red')
    axis equal 
    title(['Projeção das variaveis nos 3 principais componentes (' sprintf('%0.2g',PC3_var*100) '%)'])
    xlabel('PC1') 
    ylabel('PC2') 
    zlabel('PC3')
    legend('Primeiro Andar','Segundo Andar')
    saveas(gcf,exp_folder+'/PC3.png')
    
    figure
    scatter(ux_p(:,1),ux_p(:,2),5,'green')
    hold on;
    scatter(ux_s(:,1),ux_s(:,2),5,'red')
    axis equal
    title(['Projeção das variaveis nos 2 principais componentes (' sprintf('%0.2g',PC2_var*100) '%)'])
    xlabel('PC1') 
    ylabel('PC2')
    legend('Primeiro Andar','Segundo Andar')
    saveas(gcf,exp_folder+'/PC2.png')
    
    
    
    %disp(norm(V(:, 1)));
    %autofaces primeiro andar
    Wp = zeros(f, 1);
    
    for i = 1:PA_ex
      Wp = Wp + (V' * Phi(:, i)); % (11) -> Phi = x - u 
    end
    
    Wp = Wp / PA_ex;
    
    %autofaces segundo andar
    Ws = zeros(f, 1);
    for i = SA_off:Tex
      Ws = Ws + (V' * Phi(:, i)); % (11)
    end
    
    Ws = Ws / SA_ex;
    
    %------------------------------------TESTES------------------------------------
    % se o treino foi realizado no set initial, utiliza o set de dados
    % aumentados
    if reload == 1 || tc_row == 1
        
        if set == "initial"
           tst_set = 'augmented'
           PAt_files = get_files(['audio/' tst_set '/primeiro_andar']);
           SAt_files = get_files(['audio/' tst_set '/segundo_andar']);
           PAt_ex = length(PAt_files);
           SAt_ex = length(SAt_files);
        end
        
        if guess_other == 1
            O_files = get_files(['audio/initial/outro']);
            O_files = [O_files;get_files(['audio/augmented/outro'])];
            lbls = ["Primeiro Andar","Segundo Andar","Outro"];
            if confusion == 1
                O_files = [O_files;get_files(['audio/augmented/confusion'])];
            end
        else
            O_files = cell(0,1);
            lbls = ["Primeiro Andar","Segundo Andar"];
        end
    end
    
    O_ex = length(O_files);
    
    true_lbls = ones(1,PAt_ex); %predições corretas
    true_lbls = [true_lbls ones(1,SAt_ex).*2];
    true_lbls = [true_lbls ones(1,O_ex).*3];
    pred_lbls = zeros(1,length(true_lbls));
    
    
    Filepath_t = [PAt_files; SAt_files; O_files];
    % Filepath_t = [O_files];
    tst_ex = length(Filepath_t)
    Xt = zeros(ma,tst_ex);
    DFFS = zeros(tst_ex,1);
    
    for i = 1 : tst_ex
        %reduz dimensão das amostras
        Filepath_t{i,1}
        xt = resize(Filepath_t{i,1},ma);
        Xt(:,i) = xt;
        %-----------------------------
    
        Phit = xt' - mu;
        Wt = V' * Phit; % (11)
        DFFS(i) = (Phit' * Phit) - sum(Wt.^2); % (13)
    
        DIFSp = norm(Wt - Wp) % (14)
        DIFSs = norm(Wt - Ws) % (14)
    
        if DFFS(i) < Dp || guess_other == 0
            if DIFSp < DIFSs
                pred_lbls(i) = 1;
                fprintf('PRIMEIRO ANDAR !!!\n')
            else
                pred_lbls(i) = 2;
                fprintf('SEGUNDO ANDAR !!!\n')
            end
        else
            pred_lbls(i) = 3;
            fprintf('OUTRO !!!\n')
        end
    
    end
    
    M_DFFS = mean(DFFS)
    S_DFFS = std(DFFS)
    precision = @(confusionMat) diag(confusionMat)./sum(confusionMat,2);
    recall = @(confusionMat) diag(confusionMat)./sum(confusionMat,1)';
    
    figure
    C = confusionmat(true_lbls,pred_lbls);
    confusionchart(C,lbls);
    title('Matriz de confusão')
    Accuracy = sum(true_lbls == pred_lbls,'all')/numel(pred_lbls)
    Precision = precision(C);
    Recall = recall(C);
    Precision(isnan(Precision))=0
    Recall(isnan(Recall))=0
    saveas(gcf,exp_folder+'/confmatrix.png')
    
    figure
    title('Distancia do espaço de categorias')
    histogram(DFFS,'Normalization','pdf')
    saveas(gcf,exp_folder+'/dffs.png')
    
    %guarda as metricas do experimento
    exp_file = folder+"/experiments.csv"
    exp_file_data = [exp_no,Accuracy,Tex,tst_ex,Tst_p,p,guess_other,confusion,Dp]
    writematrix(exp_file_data,exp_file,'WriteMode','append');
    
    writematrix(C,exp_folder+"/confmatrix.csv");
    writematrix([Precision Recall],exp_folder+"/metrics.csv");
end