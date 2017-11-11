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

New in this version (Jan. 27, 2017): Also calculates areas for each peak, and exports as CSVs called Raw.Areas, Corr.Areas, Norm.Raw.Areas, and Norm.Corr.Areas (defined similarly to the peak height files). Filenames of peak height output files have also been changed to avoid confusing them with the area files. For peaks that share the same baseline, the baseline remains unchanged and the areas are defined between each endpoint and the trough between the peaks. So to get the total area of the aliphatic region (for example), add the areas of the aliph28 and aliph29 peaks.

For details on the use of this program, including a tutorial and description of all output files, see the included file "procedure for R program.docx".

## Licensing

Copyright © 2017 Suzanne Hodgkins and Florida State University.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
