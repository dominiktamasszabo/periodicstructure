function [basisSectionBegin, basisSectionEnd] = generateBasisSections()

cellBasisSize = 1.1;
cellBasisDivs = 16;
cellBasisOverlap = 1;

metalBasisSize = 0.15;
metalBasisDivs = 20;
metalBasisOverlap = 1;

[cellSquareBegin, cellSquareEnd] = generateSquareSections(cellBasisSize, cellBasisDivs, cellBasisOverlap);
[metalSquareBegin, metalSquareEnd] = generateSquareSections(metalBasisSize, metalBasisDivs, metalBasisOverlap);

basisSectionBegin = [cellSquareBegin, metalSquareBegin];
basisSectionEnd = [cellSquareEnd, metalSquareEnd];

end