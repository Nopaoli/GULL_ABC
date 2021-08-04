//Parameters for the coalescence simulation program : fastsimcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NPOP1
NPOP2
//Haploid samples sizes 
15
26
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0	MIG21	  	
MIG12  	0
//Migration matrix 1
0       0
0	0
//Migration matrix 2
0	0   
0       0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
T2 0 0 1 1 0 1
T1 0 1 1 RSANC 0 2
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
DNA 1000 0.00000 MUTRATE 0.66
