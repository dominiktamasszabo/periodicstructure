function [basisSectionBegin, basisSectionEnd] = generateBasisSections()

cellBasisSize = 1.1;
cellBasisDivs = 4;
cellBasisOverlap = 0;

metalBasisSize = 0.15;
metalBasisDivs = 1;
metalBasisOverlap = 0;

[cellSquareBegin, cellSquareEnd] = generateSquareSections(cellBasisSize, cellBasisDivs, cellBasisOverlap);
[metalSquareBegin, metalSquareEnd] = generateSquareSections(metalBasisSize, metalBasisDivs, metalBasisOverlap);

basisSectionBegin = [cellSquareBegin, metalSquareBegin];
basisSectionEnd = [cellSquareEnd, metalSquareEnd];

end