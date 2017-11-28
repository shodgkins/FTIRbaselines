# FTIRbaselines
Processes FTIR spectra of natural organic matter by finding exact peak locations, baseline-correcting the peaks, and converting them to relative abundances.

## Overview

This script finds the exact wavenumber locations of specific peaks in FTIR spectra of natural organic matter (such as plants and peat soils). It then baseline-corrects the peaks and converts them into relative abundances (relative to the integrated area of the whole spectrum). It does this for the following peaks (approx. wavenumbers in cm^-1):

carb (1030): carbohydrates
arom15 (1510): aromatics
arom16 (1630) : aromatics, or deprotonated COO-
acids (1720): protonated COOH
aliph28 (2850): aliphatics, lipids, and waxes
aliph29 (2920): aliphatics, lipids, and waxes

For details on the use of this program, including a tutorial and description of all output files, see the included file "procedure for R program.docx".

An example dataset is included in the file "example_data.csv". Humification indices for these spectra (based on fixed wavenumbers, not using this script) are published in the following manuscript:
Hodgkins S. B., Tfaily M. M., McCalley C. K., Logan T. A., Crill P. M., Saleska S. R., Rich V. I. and Chanton J. P. (2014) Changes in peat chemistry associated with permafrost thaw increase greenhouse gas production. Proc. Natl. Acad. Sci. USA 111, 5819–5824.

## Licensing

Copyright © 2017 Suzanne Hodgkins and Florida State University.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
