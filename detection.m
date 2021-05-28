close all, clear all

gt = xml2struct('PETS2009-S2L1.xml');

path = 'View_001/';

nFrames = size(dir([path '/*.jpg']), 1); %dir da-nos uma lista com os nomes de todas as frames .jpg e o size da-nos o numero de frames

% as proximas duas linhas foram calculadas aqui para conseguirmos criar o Bkg com o tamanho das frames
strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', 0, 'jpg'); %string com o nome da primeira frame
I = imread(strFrame);
Bkg = zeros(size(I)); %imagem a zeros (inicialmente) do tamanho das nossas frames

alfa = 0.005;

% for i=0:(nFrames-1) %vamos percorrer todas as frames
%     strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', i, 'jpg');
%     Y = imread(strFrame);
%     Bkg = alfa * double(Y) + (1-alfa) * double(Bkg);
%     imagesc(uint8(Bkg)); axis ij, drawnow
% end
% 
% imwrite(uint8(Bkg), 'bkg.png');

imgbk = imread('bkg.png');
imgMap = zeros(size(imgbk));

thr = 75;
minArea = 200;

se = strel('disk', 3);

figure;
for i=0:(nFrames-1) % ler frames sequencialmente e para cada imagem calcular a diferença com a imagem de background 
    strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', i, 'jpg');
    imgfr = imread(strFrame); %para ir buscar cada imagem
%     ax(1) = subplot(1,2,1); 
    imshow(imgfr); title('Pedestrian Detection'); hold on;
    
    imgdif = (abs(double(imgbk(:,:,1))-double(imgfr(:,:,1))) > thr) | (abs(double(imgbk(:,:,2))-double(imgfr(:,:,2))) > thr) | (abs(double(imgbk(:,:,3))-double(imgfr(:,:,3))) > thr);
    % imgdif só fica ativo (a 1) no sítio das onde há movimento aka onde há
    % diferenças

    bw = imclose(imgdif, se); %aplicar operação morfológica para limpara as regiões um bocadinho melhor
    % mesmo assim acusa muitas regiões ativas que são ruído - regionprops:
    
    [lb num] = bwlabel(bw);
    regionProps = regionprops(lb, 'Area', 'BoundingBox', 'FilledImage', 'Centroid');
    inds = find([regionProps.Area] > minArea); % guarda os indices das regiões que satisfazem a condição
    regnum = length(inds); % numero de regiões que nos interessam
    
    if regnum
        for j=1:regnum
            [lin col] = find(lb == inds(j)); % devolve todas as posições [y x] da região 
            upLPoint = min([lin col]); % devolve y, x
            dWindow = max([lin col]) - upLPoint + 1; % devolve height, width
            
            rectangle('Position', [fliplr(upLPoint) fliplr(dWindow)], 'EdgeColor', [1 1 0], 'linewidth', 2); %fliplr porque precisamos que position = [x y w h]
            text(regionProps(inds(j)).Centroid(1), regionProps(inds(j)).Centroid(2), num2str(j), 'Color', 'w','FontSize', 20);
        end
    end
    
    % -------------- regioes do ground truth para esta frame -------------- %
    frame = gt.Children((2*i)+2);
    %frame
    regionsList = frame.Children(2); %regionsList contem as regioes do ground truth
    %regionsList
    
    %Vamos iterar sobre todas as regiões de regionsList para pintarmos as suas boxes
    for j=2:2:length(regionsList.Children)
        region = regionsList.Children(j);
        boundingBox = region.Children(2).Attributes; % boundingBox tem o nome e o valor das coordenadas da box
        
        w = str2double(boundingBox(2).Value);
        h = str2double(boundingBox(1).Value);
        x = str2double(boundingBox(3).Value)-(w/2);
        y = str2double(boundingBox(4).Value)-(h/2);
        
        rectangle('Position', [x y w h], 'EdgeColor', [1 1 1], 'linewidth', 2);
    
    end    
    
    % --------------------------------------------------------------------- %
            
    drawnow
       
%     for k=1:length(inds) 
%             %centroid = regionProps(inds(k)).Centroid;
%             %imgMap(round(centroid(2)), round(centroid(1))) = imgMap(round(centroid(2)), round(centroid(1))) + 1;
%         [lin col] = find(lb == inds(k)); % devolve todas as posições [y x] da região
%         for pos=1:length([lin col])
%             imgMap(lin(pos), col(pos)) = imgMap(lin(pos), col(pos)) + 1;
%         end
%         
%     end
%     if rem(i,10) == 0
%         v_min = min(imgMap(:));
%         v_max = max(imgMap(:));
% 
%         ax(2) = subplot(1,2,2); imshow(imgMap); title('Heatmap'); hold on
%         colormap(ax(2), flipud(jet)); colorbar; hold on
%         caxis([v_min v_max]);
%         drawnow
%     end
    
end

