### Search for peptide synthesis side products of a certain mass

Synthesis of long peptides via solid phase peptide synthesis can often lead to a number of side products that include truncated peptides and peptides with deletions that decrease the yield of the synthesis and can be difficult to purify. This simple C++ utility is meant to help identify coupling steps that are problematic by determining which side products might be appearing in the LC-MS of the crude cleaved peptide mixture. The utility implements a simple recursive solution of the [subset sum problem](https://en.wikipedia.org/wiki/Subset_sum_problem) to find the candidate peptides.

The inputs to the program are a peptide string and the assumed molecular weight of the side product derived from LC-MS output. The program outputs all possible side products of peptide synthesis within a mass of +/- 1 of the input mass. When calculating peptide weights, we assume a C- terminal amide (mass difference -0.99) and a N-terminal acetyl group (mass difference +42.01). The tolerance or cap mass can be changed by modifying the first few lines of the function `findCandidatePeptides`:

```
double cap_mass = 41.0265;
double tolerance = 1.0;
```

#### Installation and Usage

Compilation:
In the command line:
```
g++ -std=c++11 peptidechecker.cpp
```

Execution:

```
./a.out
```

Example output:

```
Enter the peptide string:
```
