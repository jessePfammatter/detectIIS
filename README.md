# IISDetector

Software for the detection of interictal spikes within mouse EEG records. These scripts are associated with the paper 'A Predictive Epilepsy Index Based on Probablistic Classification of Interictal Spike Waveforms' by Jesse A. Pfammatter et al. 2018. Sample data (and the data associated with the manuscript) can be found at https://www.ebi.ac.uk/biostudies/studies/S-BSST180. These data are already in matlab stuct format -- original EDF files available upon request (jesse.pfammatter@gmail.com). 

This set of software is currently more developer friendly than end user friendly. The develop_IIS_detection.m script is the main script that calls other functions within this repository -- using this script and follow along with the manuscript in order to get things up and running.

# Associated Scripts

## build_dataFrame_falconHawk_detectionIIS_512.m

A script that was used to build that dataFrame available on the biostudies database linked above. This script isn't really useful to those that don't have the original EDF files, but I've included it as a record of how we downsampled/upsampled, filtered, and otherwise preprocessed our data.

## detectEvents_IIS.m

A shell script that calls twoThresholdPeakDetect to identify high amplitude events and packages them together into a nice format that can be used downstream.

## develop_IIS_detection.m

This is the main script that can produce a useful output for analysis etc. It's the script we used to produce all the figures for the publication.

## highPassChebyshev1Filt_EEG.m

A script that performs high pass filtering on EEG data.

## normalizeEEG.m

This script normalizes the EEG using a modification on the mean and standard deviation normalization. Rather than using the mean and standard deviation of the data we use the mean and standard of a Gaussian model fit to the data. For details please see Pfammatter et al 2018.

## read_EDF.m

A script used to read EDF files. Again, you won't need this if you start with data from the BioStudies repository linked above.

## sixtyHzFilt_EEG.m

A notch filter to remove 60Hz line noise.

## twoThresholdPeakDetect.m

A script that identifies peaks using a two threshold method. A peak is identified when the EEG signal first crosses a high threshold followed by a low threshold.



