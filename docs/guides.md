# Guide Sequences
The workflow requires a guide sequence reference to run. For each sample sequenced, the samples.csv should encode a tab separated text file that contains the guide sequences to be detected per cell (the guides column references this file by name). When launching the workflow, a folder containing each of these files referenced is passed by path to the guidesDir paramter in the workflow or [runParams.yml](examples/runParams.yml) 

Each `guides.txt` file includes two columns in tab-separated format, without column headers see [guides.txt](examples/guides.txt).

* the first column values should contain the A,T,G,C sequence of the guide to be detected
* the second column values should containg the desired name corresponding to each sequence, which will be used for the final UMI and Read matrix generation.

## Guide Structure Requirements
These guide sequences need to be of the same length, and the pipeline will enforce this by trimming them to the lowest length of the sequences provided if they are not all identical. It is important to keep this in mind as if in the process of the trimming out basepairs the sequences are no longer unique, the pipeline will still proceed. 

In addition, the pipeline requires that the guide sequences to be detected are found in the same starting position within Read 2 for the analysis meaning that variability among the constructs that would induce bp shifts in the final location will cause potential issues. 