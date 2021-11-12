function [] = PrintResultsAsFigure(TxTCell)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


FigH = figure('menubar', 'none');
ListH = uicontrol('Style', 'listbox', ...
   'Units',    'normalized', ...
   'Position', [0,0,1,1], ...
   'String',   {}, ...
   'Min', 0, 'Max', 2, ...
   'Value', [],...
   'FontSize',12);
for i = 1:length(TxTCell)
    pause(0.1);
    newString = cat(1, get(ListH, 'String'), TxTCell(i));
    set(ListH, 'String', newString);
    drawnow;
end


end

