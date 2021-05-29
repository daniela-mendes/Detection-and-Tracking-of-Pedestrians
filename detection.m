close all, clear all

%gt = xml2struct('PETS2009-S2L1.xml');
last_fr = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]}; % este array vai conter todos os centroides das regioes das ultimas 15 frames - usado para traçar trajetorias dinamicas

% ------ arrays para calcular labels ------ %
bbox_last_fr = {}; % este array vai conter as bounding boxes de todas as regioes da ultima frame visitada
idx_last_fr = []; % array com indices das regioes da ultima frame visitada
bbox_curr_fr = {}; % este array vai conter as bounding boxes de todas as regioes da frame em que estamos
idx_curr_fr = []; % este array vai conter os indices das regioes da bbox_last_fr com as quais cada umas das novas regioes interseta
% --------------------------------------- %

path = 'View_001/';

nFrames = size(dir([path '/*.jpg']), 1); %dir da-nos uma lista com os nomes de todas as frames .jpg e o size da-nos o numero de frames
I_over_O = zeros(nFrames, 1); % para guardar as intersections over unions das regions de todas as frames

% as proximas duas linhas foram calculadas aqui para conseguirmos criar o Bkg com o tamanho das frames
strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', 0, 'jpg'); %string com o nome da primeira frame
I = imread(strFrame);
Bkg = zeros(size(I)); %imagem a zeros (inicialmente) do tamanho das nossas frames

% ---------------------- Calculo da background image ---------------------- %
alfa = 0.005;

% for i=0:(nFrames-1) %vamos percorrer todas as frames
%     strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', i, 'jpg');
%     Y = imread(strFrame);
%     Bkg = alfa * double(Y) + (1-alfa) * double(Bkg);
%     imagesc(uint8(Bkg)); axis ij, drawnow
% end
% 
% imwrite(uint8(Bkg), 'bkg.png');
% ------------------------------------------------------------------------- %

imgbk = imread('bkg.png');
imgMap = zeros(size(imgbk, 1), size(imgbk, 2));
imgTraj = imread('bkg.png');

thr = 75;
minArea = 200;

se = strel('disk', 3);

figure;
for i=0:(nFrames-1) % ler frames sequencialmente e para cada imagem calcular a diferença com a imagem de background 
    
    strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', i, 'jpg');
    imgfr = imread(strFrame); %para ir buscar cada imagem
    subplot(1,2,1); imshow(imgfr); title('Pedestrian Detection'); hold on;
    
    
    % -------------- regioes do ground truth para esta frame -------------- %
%     gt_regs = zeros(20,4); % vetor onde guardamos as bounding boxes de cada regiao do ground truth
%     n=0;
%     
%     frame = gt.Children((2*i)+2);
%     regionsList = frame.Children(2); %regionsList contem as regioes do ground truth
%     
%     %Vamos iterar sobre todas as regiões de regionsList para pintarmos as suas boxes
%     for j=2:2:length(regionsList.Children)
%         region = regionsList.Children(j);
%         boundingBox = region.Children(2).Attributes; % boundingBox tem o nome e o valor das coordenadas da box
%         
%         w = str2double(boundingBox(2).Value);
%         h = str2double(boundingBox(1).Value);
%         x = str2double(boundingBox(3).Value)-(w/2);
%         y = str2double(boundingBox(4).Value)-(h/2);
%         
%         rectangle('Position', [x y w h], 'EdgeColor', [1 1 0], 'linewidth', 2);        
%         
%         gt_regs(n+1, 1) = x;
%         gt_regs(n+1, 2) = y;
%         gt_regs(n+1, 3) = w;
%         gt_regs(n+1, 4) = h;
%         n = n+1;        
%     
%     end    
% 
%     gt_regs = gt_regs(1:n,:);
    % --------------------------------------------------------------------- %
    
    
    % ------------- regioes detetadas por nós para esta frame ------------- %
    imgdif = (abs(double(imgbk(:,:,1))-double(imgfr(:,:,1))) > thr) | (abs(double(imgbk(:,:,2))-double(imgfr(:,:,2))) > thr) | (abs(double(imgbk(:,:,3))-double(imgfr(:,:,3))) > thr);
    % imgdif só fica ativo (a 1) no sítio das onde há movimento aka onde há
    % diferenças

    bw = imclose(imgdif, se); %aplicar operação morfológica para limpara as regiões um bocadinho melhor
    % mesmo assim acusa muitas regiões ativas que são ruído - regionprops:
    
    [lb num] = bwlabel(bw);
    regionProps = regionprops(lb, 'Area', 'BoundingBox', 'FilledImage', 'Centroid');
    inds = find([regionProps.Area] > minArea); % guarda os indices das regiões que satisfazem a condição
    
    last_fr{rem(i,15)+1} = []; % rem(i,15)+1=[] eh a forma de acedermos ciclicamente ah posicao do array que foi atualizada ha mais tempo e limpar os centroids associados a essa frame; 
    
    bbox_curr_fr = []; %vamos guardar as bounding boxes das regioes desta frame
    idx_curr_fr = []; %e aqui guardamos os indices com as quais essas bounding boxes intersetam
    
    for j=1:length(inds)
        [lin col] = find(lb == inds(j)); % devolve todas as posições [y x] da região 
        upLPoint = min([lin col]); % devolve y, x
        dWindow = max([lin col]) - upLPoint + 1; % devolve height, width
        box = [fliplr(upLPoint) fliplr(dWindow)];
        
        if i ~= 0 % se estivermos na primeira frame, vamos simplesmente guardar as regioes, nao temos que compara-las com as da frame anterior, porque nao existe frame anterior
            overlap = false;
            maximo = 0;
            idx_maximo = 0;
            idx_r = 0;
            for r=1:length(bbox_last_fr)
                if rectint(box, bbox_last_fr{r}) > maximo %interseccao ~= o quer dizer que ha overlap das regioes e que eh a mesmas regiao ativa
                    overlap = true;
                    maximo = rectint(box, bbox_last_fr{r});
                    idx_maximo = idx_last_fr(r); %quando ha intersecao, indices deverao ser iguais pois terao a mesma label
                    idx_r = r; % queremos guardar o r para, no calculo da trajetoria, conseguirmos ir buscar a regiao em bbox_last_fr que esta na posicao r
                end
            end
            
            if overlap == false 
                if j ~= 1 % a nossa regiao detetada nao intersetou com nenhuma regiao da frame anterior, o que significa que eh uma regiao nova
                    idx_maximo = max(max(idx_last_fr), max(idx_curr_fr)) + 1; % a nova regiao fica com uma nova label, imediatamente a seguir ah label mais alta que ja existia
                else
                    idx_maximo = max(idx_last_fr) + 1;
                end
%             else
%                 last = bbox_last_fr{idx_r}; % last eh a regiao na bbox_last_fr que eh analoga ah regiao em questao (i.e. box)
%                 centroid_x = [(last(1)+(last(3)/2)) (box(1)+(box(3)/2))];
%                 centroid_y = [(last(2)+(last(4)/2)) (box(2)+(box(4)/2))];
%                 plot(centroid_x, centroid_y, 'w-');
            end
            
            bbox_curr_fr{end+1} = box;
            idx_curr_fr(end+1) = idx_maximo;
            text(regionProps(inds(j)).Centroid(1), regionProps(inds(j)).Centroid(2)-(regionProps(inds(j)).BoundingBox(4)/2)-10, num2str(idx_curr_fr(end)), 'Color', [0.949 0.949 0.949],'FontSize', 20);
            
        else
            bbox_curr_fr{end+1} = box; %para a primeira frame, apenas guardamos bounding box de cada regiao
            idx_curr_fr(end+1) = j; %labels na primeira frame correspondem ah ordem em que vemos as regioes
            text(regionProps(inds(j)).Centroid(1), regionProps(inds(j)).Centroid(2)-(regionProps(inds(j)).BoundingBox(4)/2)-10, num2str(j), 'Color', [0.949 0.949 0.949],'FontSize', 20);
        end

        rectangle('Position', box, 'EdgeColor', [1 1 1], 'linewidth', 2); %fliplr porque precisamos que position = [x y w h]
        
        last_fr{rem(i,15)+1}(end+1) = regionProps(inds(j)).Centroid(1); % onde limpamos os centroides da frame mais antiga, escrevemos agora os centroides da frame mais recente
        last_fr{rem(i,15)+1}(end+1) = regionProps(inds(j)).Centroid(2);
        
        % ----------- IoU ----------- %
%         for d=1:size(gt_regs, 1) % gt_regs tem as bounding boxes das regions do ground truth e queremos iterar sobre cada regiao do gt para calcular a intersecao com a regiao detetada por nós
%             box_gt = [gt_regs(d, 1) gt_regs(d, 2) gt_regs(d, 3) gt_regs(d, 4)];
%             intersection = rectint(box, box_gt); %rectint calcula a area de intersecao entre as duas regioes dadas
%             union = (box(3) * box(4)) + (box_gt(3) * box_gt(4)) - intersection;
%             i_o_u = intersection/union;
%             % if i_o_U ~=0; I_over_U !!!!!! bboxOverlapRatio
%         end
        % --------------------------- %
        
    end
    
    bbox_last_fr = bbox_curr_fr; %
    idx_last_fr = idx_curr_fr;
    
    n_fr = min(i+1, 15); %para evitar que, nas primeiras 3 frames, tentemos acessar os centroides de frames que ainda não visitamos
    for f=1:n_fr
        plot(last_fr{f}([1:2:length(last_fr{f})]), last_fr{f}([2:2:length(last_fr{f})]), 'w*'); % plot das trajetorias dinâmicas
    end
    
    drawnow
    % --------------------------------------------------------------------- %
    
    
    % ------------- trajetorias realizadas pelos pedestres ---------------- %
    for k=1:length(inds) 
        centroid_x = round(regionProps(inds(k)).Centroid(1));
        centroid_y = round(regionProps(inds(k)).Centroid(2));
        imgTraj(centroid_y-1:centroid_y+1, centroid_x-1:centroid_x+1, 1) = 255;
        imgTraj(centroid_y-1:centroid_y+1, centroid_x-1:centroid_x+1, 2:3) = 0;
    end
    subplot(1,2,2); imshow(imgTraj); title('Performed Trajectories');
    % --------------------------------------------------------------------- %
    
end