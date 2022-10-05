# Optogenetic frequency scrambling of hippocampal theta oscillations dissociates working memory retrieval from hippocampal spatiotemporal codes
Main codebase for running analyses. Link to manuscript: https://doi.org/10.1101/2021.12.24.474139

## System requirements
- Matlab > 2021b
This has been tested on MacOS11 and 12 but should work on other Unix systems (Linux) and Windows with no changes

## Installation
Clone this repository (or download a zip file - not recommended for version control). Not further installation is required

## Usage
Each folder (behavior, CalciumImaging, and Electrophysiology) contains dedicated scripts.
Matlab files (.m) are used to extract calcium imaging data as well as perform spectral analyses (Fourier and Wavelet transforms).

For most calcium imaging functions, parameters can be set in the script header, and the script can simply be executed while being in the working folder.
Example:
```
batchDecode_PC_LMTdata
```
to decode spatial information for all mice (no arguments needed). Expected batch functions runtime ~5h on a M1 mac. Other functions run under a minute.

Functions related to electrophysiology take 'data' as input, which is corresponds to the singular datastructure (Matlab cell) containing all informations (channel names, TTL timestamps, ...) related to a given recording. In contrast to calcium imaging functions, these scripts simply require a reference to the data structure as the only argument.
Example:
```
results = shortStimsAnalysis(data)
```
where data is the data structure containing all electrophysiological recordings for one session, and results will contain all the results of the analysis.

## Contribute
For small changes to the codebase, feel free to submit a pull request. For bigger changes or issues, submit an issue.

## License
This project is licensed under the terms of the MIT license:

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


