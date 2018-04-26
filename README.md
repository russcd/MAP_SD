# MAP_SD

This software takes paired allele counts from somatic and germline samples to estimate the distortion 
effect size and the map position of segregation distortion elements in the somtic sample. 

It takes two input files.
1. A data file with read counts supporting either parental chromosome by site. 
2. A file that specifies a single chromsome-wide recombination rate in Morgans/bp for each chromosome to be considered. 

The first file is tab delimited, with no header line and has the format:
1. Chromosome ID. 
2. Position in basepairs.
3. Map position in centimorgans.
3. Read counts for the first parental type in the somatic library. 
4. Read counts for the second parental type in the somatic library. 
5. Read counts for the first parental type in the gamete library. 
6. Read counts for the second parental type in the gamete library. 

Basic Usage:

map_sd -d <data_file> 

This will produce an output file with the following format
1. Chromosome ID
2. Distorter Position
3. Ancestry ratio in the somatic library. 
4. Ancestry ratio in the gamete library. 
5. lnL of somatic ratio in gamete library. 
6. lnL of the gamete ratio in the gamete library. 

Additional Options/Advanced Usage

-w [int]    will produce windowed estimates of SD across each chromosome with window size -w <int>
  
-b [int]    will produce bootstrap estimates of SD position via resampling
  
-e [float]  uniform error rate per read
  
-t [float]  terminate optimation if successive points are within <float> likelihood units. 
