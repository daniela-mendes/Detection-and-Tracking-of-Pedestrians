close all, clear all

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
for i=0:(nFrames-1) % ler frames sequencialmente e para cada imagem calcular a diferen�a com a imagem de background 
    strFrame = sprintf('%s%s%.4d.%s', path, 'frame_', i, 'jpg');
    imgfr = imread(strFrame); %para ir buscar cada imagem
    subplot(1,2,1); imshow(imgfr); title('Pedestrian Detection'); hold on;
    
    imgdif = (abs(double(imgbk(:,:,1))-double(imgfr(:,:,1))) > thr) | (abs(double(imgbk(:,:,2))-double(imgfr(:,:,2))) > thr) | (abs(double(imgbk(:,:,3))-double(imgfr(:,:,3))) > thr);
    % imgdif s� fica ativo (a 1) no s�tio das onde h� movimento aka onde h�
    % diferen�as

    bw = imclose(imgdif, se); %aplicar opera��o morfol�gica para limpara as regi�es um bocadinho melhor
    % mesmo assim acusa muitas regi�es ativas que s�o ru�do - regionprops:
    
    [lb num] = bwlabel(bw);
    regionProps = regionprops(lb, 'Area', 'BoundingBox', 'FilledImage', 'Centroid');
    inds = find([regionProps.Area] > minArea); % guarda os indices das regi�es que satisfazem a condi��o
    regnum = length(inds); % numero de regi�es que nos interessam
    
    if regnum
        for j=1:regnum
            [lin col] = find(lb == inds(j));
            upLPoint = min([lin col]);
            dWindow = max([lin col]) - upLPoint + 1;
            
            rectangle('Position', [fliplr(upLPoint) fliplr(dWindow)], 'EdgeColor', [1 1 0], 'linewidth', 2);
            text(regionProps(inds(j)).Centroid(1), regionProps(inds(j)).Centroid(2), num2str(j), 'Color', 'w','FontSize', 20);
        end
    end
            
    drawnow
       
    for k=1:length(inds)
        centroid = regionProps(inds(k)).Centroid;
        
        for l=1:size(imgMap,1)
            for c=1:size(imgMap,2)
                dist = sqrt((c - centroid(1)).^2 + (l - centroid(2)).^2);
                if (k == 1 && i == 0)
                    imgMap(l, c) = dist;
                else
                    if imgMap(l, c) > dist
                        imgMap(l, c) = dist;
                    end
                end
            end
        end
        
    end
        
    if i == 0
        imgMap
    end
    subplot(1,2,2); imshow(imgMap); title('Heatmap');
    
end

