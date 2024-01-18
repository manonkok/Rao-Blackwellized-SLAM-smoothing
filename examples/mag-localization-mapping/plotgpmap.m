function [clims]=plotgpmap(dEft,dVarft,G1,G2,Y1,Y2,path,y,filename,component)

  % Background
  %filename = 'frame.png';

  % Figure
  hold on
  
  % Load image
  if ~isempty(filename)
    img = imread(filename);
    imagesc(img)
  end

  % Show grid
  plot(G1,G2,'-',G1',G2','-','Color',[.7 .7 .7],'LineWidth',.5)
 
  % Reshape magnetic map
  if numel(component)==1
    map = reshape(dEft(:,component),size(Y1));  
  else
    map = reshape(sqrt(sum(dEft(:,component).^2,2)),size(Y1));
  end
  map(isnan(map)) = 0;    
  
  % Reshape alpha channel
  alphamap = exp(-reshape(sqrt(max(dVarft(:,component),[],2)),size(Y1)));    
  alphamap(isnan(alphamap)) = 0;

  % Show surface
  h=pcolor(Y1,Y2,map);
  set(h,'EdgeColor','none','FaceColor','interp', ...
    'AlphaData',alphamap,'AlphaDataMapping','scaled','FaceAlpha','interp');
  
  % Show vehicle path
  plot(path(:,1),path(:,2),'--k')
  
  % Show measurements
  if numel(component)==1
    foo = y(:,component);
  else
    foo = sqrt(sum(y(:,component).^2,2));
  end
  clims = [min(foo(:)),max(foo(:))];
  %scatter(p(:,1),p(:,2),12,foo,'filled')
  
  % Color axis
  caxis(clims)
  colormap(parula)

  % Axis
  axis ij image off
  xlim([0 1920]), ylim([0 1080])

  % BG Color
  set(gcf,'Color','w')
  