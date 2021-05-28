gt = xml2struct('PETS2009-S2L1.xml');
%gt = dataset, gt.Children = frames, gt.Children.Attributes = nome do
%atributo de cada frame (number) e value (que é o number da frame)

gt.Children(2).Children(2).Children(2).Children(2).Attributes


% gt.Children(2).Children(2).Children(2).Children(2) - 1a box (region)
% gt.Children(2).Children(2).Children(2) - 1o object (1 object para cada box)
% gt.Children(2).Children(2) - 1a object list (tem todos os objects da
% frame)
% gt.Children(i+2) - 1a frame do dataset
