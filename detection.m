close all, clear all

gt = xml2struct('PETS2009-S2L1.xml');
last_fr = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]}; % este array vai conter todos os centroides das regioes das ultimas 15 frames

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

thr = 75;
minArea = 200;

se = strel('disk', 3);

figure;
for i=0:(nFrames-1) % ler frames sequencialmente e para cada imagem calcular a diferen�a com a imagem de background 
    
    strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', i, 'jpg');
    imgfr = imread(strFrame); %para ir buscar cada imagem
    ax(1) = subplot(1,2,1); imshow(imgfr); title('Pedestrian Detection'); hold on;
    
    
    % -------------- regioes do ground truth para esta frame -------------- %
    gt_regs = zeros(20,4); % vetor onde guardamos as bounding boxes de cada regiao do ground truth
    n=0;
    
    frame = gt.Children((2*i)+2);
    regionsList = frame.Children(2); %regionsList contem as regioes do ground truth
    
    %Vamos iterar sobre todas as regi�es de regionsList para pintarmos as suas boxes
    for j=2:2:length(regionsList.Children)
        region = regionsList.Children(j);
        boundingBox = region.Children(2).Attributes; % boundingBox tem o nome e o valor das coordenadas da box
        
        w = str2double(boundingBox(2).Value);
        h = str2double(boundingBox(1).Value);
        x = str2double(boundingBox(3).Value)-(w/2);
        y = str2double(boundingBox(4).Value)-(h/2);
        
        rectangle('Position', [x y w h], 'EdgeColor', [1 1 1], 'linewidth', 2);        
        
        gt_regs(n+1, 1) = x;
        gt_regs(n+1, 2) = y;
        gt_regs(n+1, 3) = w;
        gt_regs(n+1, 4) = h;
        n = n+1;        
    
    end    

    gt_regs = gt_regs(1:n,:);
    % --------------------------------------------------------------------- %
    
    
    % ------------- regioes detetadas por n�s para esta frame ------------- %
    imgdif = (abs(double(imgbk(:,:,1))-double(imgfr(:,:,1))) > thr) | (abs(double(imgbk(:,:,2))-double(imgfr(:,:,2))) > thr) | (abs(double(imgbk(:,:,3))-double(imgfr(:,:,3))) > thr);
    % imgdif s� fica ativo (a 1) no s�tio das onde h� movimento aka onde h�
    % diferen�as

    bw = imclose(imgdif, se); %aplicar opera��o morfol�gica para limpara as regi�es um bocadinho melhor
    % mesmo assim acusa muitas regi�es ativas que s�o ru�do - regionprops:
    
    [lb num] = bwlabel(bw);
    regionProps = regionprops(lb, 'Area', 'BoundingBox', 'FilledImage', 'Centroid');
    inds = find([regionProps.Area] > minArea); % guarda os indices das regi�es que satisfazem a condi��o
    
    last_fr{rem(i,15)+1} = []; % rem(i,15)+1=[] eh a forma de acedermos ciclicamente ah posicao do array que foi atualizada ha mais tempo e limpar os centroids associados a essa frame; 
    for j=1:length(inds)
        [lin col] = find(lb == inds(j)); % devolve todas as posi��es [y x] da regi�o 
        upLPoint = min([lin col]); % devolve y, x
        dWindow = max([lin col]) - upLPoint + 1; % devolve height, width
        box = [fliplr(upLPoint) fliplr(dWindow)];

        rectangle('Position', box, 'EdgeColor', [1 1 0], 'linewidth', 2); %fliplr porque precisamos que position = [x y w h]
        text(regionProps(inds(j)).Centroid(1), regionProps(inds(j)).Centroid(2), num2str(j), 'Color', 'w','FontSize', 20);
        
        last_fr{rem(i,15)+1}(end+1) = regionProps(inds(j)).Centroid(1); % onde limpamos os centroides da frame mais antiga, escrevemos agora os centroides da frame mais recente
        last_fr{rem(i,15)+1}(end+1) = regionProps(inds(j)).Centroid(2);
        
        % ----------- IoU ----------- %
        for d=1:size(gt_regs, 1) % gt_regs tem as bounding boxes das regions do ground truth e queremos iterar sobre cada regiao do gt para calcular a intersecao com a regiao detetada por n�s
            box_gt = [gt_regs(d, 1) gt_regs(d, 2) gt_regs(d, 3) gt_regs(d, 4)];
            intersection = rectint(box, box_gt); %rectint calcula a area de intersecao entre as duas regioes dadas
            union = (box(3) * box(4)) + (box_gt(3) * box_gt(4)) - intersection;
            i_o_u = intersection/union;
            % if i_o_U ~=0; I_over_U !!!!!! bboxOverlapRatio
        end
        % --------------------------- %
        
    end
    
    n_fr = min(i+1, 15); %para evitar que, nas primeiras 3 frames, tentemos acessar os centroides de frames que ainda n�o visitamos
    for f=1:n_fr
        plot(last_fr{f}([1:2:length(last_fr{f})]), last_fr{f}([2:2:length(last_fr{f})]), 'y*'); % plot das trajetorias
    end
    
    drawnow
    % --------------------------------------------------------------------- %
    
    
    % --------------------------- heatmap --------------------------------- %   
    for k=1:length(inds) 
            %centroid = regionProps(inds(k)).Centroid;
            %imgMap(round(centroid(2)), round(centroid(1))) = imgMap(round(centroid(2)), round(centroid(1))) + 1;
        [lin col] = find(lb == inds(k)); % devolve todas as posi��es [y x] da regi�o
        for pos=1:length([lin col])
            imgMap(lin(pos), col(pos)) = imgMap(lin(pos), col(pos)) + 1;
        end
        
    end
    if rem(i,10) == 0
        v_min = min(imgMap(:));
        v_max = max(imgMap(:));

        ax(2) = subplot(1,2,2); imshow(imgMap); title('Heatmap'); hold on
        colormap(ax(2), jet); colorbar; hold on
        caxis([v_min v_max]);
        drawnow
    end
    % --------------------------------------------------------------------- %
    
end

% ---------------------- pergunta do optical flow ------------------------- %
fprintf('%s\n%s\n', 'Podemos colecionar alguns centr�ides de modo a estimar um vetor que nos dar� informa��es sobre a dire��o, trajet�ria e velocidade por exemplo, tendo em conta que a footage se trata de c�mera fixa.', 'Desta forma podemos criar estimativas de onde se ir� encontrar um objeto numa dada frame, desde que a sua trajet�ria n�o sofra mudan�a de dire��o. V=centroide_atual-centroide_anterior/intervalo de tempo (de uma frame para outra por exemplo)');
% ------------------------------------------------------------------------- %

