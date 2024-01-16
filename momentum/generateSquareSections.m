function [begins, ends] = generateSquareSections(size, divs, overlap)

assert (overlap>=0, "Negativ overlap") ;

borders = linspace(-size/2, +size/2, divs + overlap + 1);

% upperCellBorder
upperBorderBeginX = borders(1:divs);
upperBorderBeginY = repmat(size/2, 1, divs);
upperBorderEndX = borders(2 + overlap:divs + overlap + 1);
upperBorderEndY = repmat(size/2, 1, divs);

upperBorderBegin = [upperBorderBeginX;upperBorderBeginY];
upperBorderEnd = [upperBorderEndX; upperBorderEndY];

% lowerBorder

lowerBorderBegin = [upperBorderBeginX;-upperBorderBeginY];
lowerBorderEnd = [upperBorderEndX; -upperBorderEndY];

% leftBorder

leftBorderBegin = [-upperBorderBeginY; upperBorderBeginX];
leftBorderEnd = [-upperBorderEndY; upperBorderEndX];

%rightBorder

rightBorderBegin = [upperBorderBeginY; upperBorderBeginX];
rightBorderEnd = [upperBorderEndY; upperBorderEndX];

% Merge
begins = [upperBorderBegin, lowerBorderBegin, leftBorderBegin, rightBorderBegin];
ends = [upperBorderEnd, lowerBorderEnd, leftBorderEnd, rightBorderEnd];


end