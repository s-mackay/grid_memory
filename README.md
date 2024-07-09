# Source data and code

This repository contains the source data and code to regenerate figures for "Concept and location neurons in the human brain provide the what and where in memory formation"

## Matlab version

This code was developed on Matlab R2019b, but should also run on later versions. Minor adjustments may be required.

## How to run

While There are separate scripts for each figure, I recommend running them from the batch script 

all_figures_and_stats.m

If you do not provide any arguments you will be prompted for the path to this git repo and whether you would like to read .mat or .xlsx files.
The easiest way to get started would be to add the files to your matlab path and then run the all_figures_and_stats
```matlab
data_directory = '/home/jdoe/grid_memory';
read_from_excel = false;
addpath(genpath(data_directory));
all_figures_and_stats(data_directory, read_from_excel);
```

## .mat and .xlsx files

As per requirement of the journal, I tried to provide all source files both in .mat and .xlsx format. This was not possible for figure 3, S3 and S4, since their source files were to large to be included in a git repo when saved as .xlsx.
Reading .xlsx files is slow and not recommended.
