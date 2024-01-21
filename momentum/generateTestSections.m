function [potTestSectionBegin, potTestSectionEnd, forceTestSectionBegin, forceTestSectionEnd, sints] = generateTestSections()

global cellSize;
global metalSize;
global voltage;

cellTestSize = cellSize;
cellTestDivs = 24;
cellTestOverlap = 0;

metalTestSize = metalSize;
metalTestDivs = 48;
metalTestOverlap = 0;

[cellSquareBegin, cellSquareEnd] = generateSquareSections(cellTestSize, cellTestDivs, cellTestOverlap);
[metalSquareBegin, metalSquareEnd] = generateSquareSections(metalTestSize, metalTestDivs, metalTestOverlap);

testSectionBegin = [cellSquareBegin, metalSquareBegin];
testSectionEnd = [cellSquareEnd, metalSquareEnd];

% Force testSections
forceTestSectionBegin = testSectionBegin(:, 1:cellTestDivs*2);
forceTestSectionEnd = testSectionEnd(:, 1:cellTestDivs*2);

% Pot testSections
potTestSectionBegin = testSectionBegin(:, cellTestDivs*2+1:length(testSectionBegin));
potTestSectionEnd = testSectionEnd(:, cellTestDivs*2+1:length(testSectionBegin));

% Integral values to solve for

potIntValues = zeros(length(potTestSectionBegin), 1);
% Left
potIntValues(1:cellTestDivs) = 2 * voltage / 2;
% Right
potIntValues(cellTestDivs+1:2*cellTestDivs) = -1 * potIntValues(1:cellTestDivs);

% Metal
potIntValues(2*cellTestDivs + 1:2*cellTestDivs+4*metalTestDivs) = zeros(4*metalTestDivs,1);
forceIntValues = zeros(2*cellTestDivs,1);

sints = [potIntValues; forceIntValues];

end