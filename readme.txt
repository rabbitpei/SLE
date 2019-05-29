MATLAB environment.

Two sample dataset:

LUAD input data:
PPI network: adj_edges_all.txt
microarray data: pruning_normal_expression_step2.txt, pruning_tumor_expression_step2.txt

FLU input data: 
PPI network: adj_edges_all.txt
gene name: genename.mat
microarray data: T1.mat

running pipeline:
firstly run: *_step1.m
secondly run: *_step2.m

permutate the n reference samples:
firstly run: *_step1_permutation.m
secondly run: *_step2_permutation.m
