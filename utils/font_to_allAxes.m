% Apply the publication's Arial font to axes and text boxes in the current figure.
allAxesInFigure = findall(gcf, 'type', 'axes');
allTextBoxes = findall(gcf, 'type', 'textbox');
for jj = 1:length(allAxesInFigure)
    set(allAxesInFigure(jj), 'fontname', 'Arial')
end

for jj = 1:length(allTextBoxes)
    set(allTextBoxes(jj), 'fontname', 'Arial')
end
