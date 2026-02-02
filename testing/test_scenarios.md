# Evaluation of Tracking Methods in Julian
Here, we employ the tracking technique implemented in Julian to demonstrate its capability and compare it with previously developed MATLAB codes.

## Key Features
Compared with the MATLAB code, the implementation based on [MeshArray](https://juliaclimate.github.io/MeshArrays.jl/dev/) enables the toolbox to handle cross-boundary problems (e.g., longitude wrapping from 0°E to 360°E), which are relatively difficult to address in the MATLAB environment. Idealized and real-world examples are provided below.

## Idealized Examples
### General Scenario
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_see.gif)
### Spliting Scenario
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_spliting.gif)
### Merging Scenario
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_merging.gif)
### Complex Scenario
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_60.gif)

## Comparing with MATLAB codes
To test if the model can provide the same outputs as what we gets from the original MATLAB code, we manually seperate the 0.25o NOAA OI SST data during 1982-2022 in the Southeast Australia (Lon-Lat-Time; 400-160-14975;data_old) into two parts.
![two_parts](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/face_example_aus.png)

This results in a new dataset with dimensions 200–160–2–14,975 (denoted as `data_new`). We then perform tracking on both `data_old` and `data_new` using MATLAB and Julia codes. All computations are carried out on the NCAR HPC system, with jobs submitted using a single CPU. Outputs can be found in githubs for MATLAB ([.mat](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/tracks_example.mat)) and Julia ([.h5](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/tracks_example.h5.zip))

<table>
<colgroup>
<col width="10%" />
<col width="45%" />
<col width="45%" />
</colgroup>
<thead>
<tr class="header">
<th></th>
<th>Julia</th>
<th>MATLAB</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Number of Events</td>
<td>984</td>
<td>984</td>
</tr>
<tr class="even">
<td>Spliting/Merging Ratio</td>
<td>9.55%</td>
<td>9.55%</td>
</tr>
<tr class="odd">
<td>Number of Objects</td>
<td>20254</td>
<td>20254</td>
</tr>
<tr class="even">
<td>Number of Pixels</td>
<td>37148175</td>
<td>37148175</td>
</tr>
<tr class="odd">
<td>Running Time</td>
<td>0:30:02</td>
<td>0:32:13</td>
</tr>
</tbody>
</table>

### The 2011 WA MHWs in MATLAB and Julia Codes - The one names MHWs
![see](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhwbranknn2.gif)

## Global Tracks in ECCOV4R4
Here we do the MHW tracking in ECCOV4R4 during 1992-2017 (90-90-13-9497; Lon-Lat-Face-Time) with MHW information identified based on baselines during 1992-2012, following the Hobday et al. (2016) definition. Key statistics are provided below:

<table>
<colgroup>
<col width="17%" />
<col width="83%" />
</colgroup>
<thead>
<tr class="header">
<th>Metric</th>
<th>Value</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Number of Events</td>
<td>5217</td>
</tr>
<tr class="even">
<td>Spliting/Merging Ratio</td>
<td>12.1%</td>
</tr>
</tbody>
</table>

Some global statistics are also provided here:
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/Figure1_ce.png)

Some historically recongized events 
### The "Blob"
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_blob_ecco_2.gif)

### The 1997-98 ENSO
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_enso_ecco_2.gif)

### A Central Pacific Nino during 1997-98
![test](https://github.com/ZijieZhaoMMHW/Tracking_Julian/blob/main/mhw_cp_ecco_2.gif)

## A Detailed Analysis for MHW tracks in Australia Regions as an example
We then do some detailed analysis on the 984 Southern Australia Examples.


