MIToS
===== 
(Mutual Information Tools for protein Sequence analysis ported from [MIToS.jl](https://github.com/diegozea/MIToS.jl))

MIToS is an environment for Mutual Information (MI) analysis and implements several useful tools for Multiple Sequence Alignments (MSAs) and PDB structures management in the Python language. MI allows determining covariation between positions in a MSA. MI derived scores are good predictors of residue contacts and functional sites in proteins [1,2].

MIToS starting point was an improvement of the algorithm published by Buslje et. al. [1]. A BLOSUM62-based pseudocount strategy, similar to Altschul et. al. [3], was implemented for a better performance in the range of MSAs with low number of sequences. **MIToS** offers all the necessary tools for using, developing and testing MI based scores, in different modules:

* **MSA** defines multiple functions and types for dealing with MSAs:
  * `AnnotatedMultipleSequenceAlignment` is a type for saving MSAs and their  `Annotations`, it is useful for working on MSAs in `Stockholm` format.
  * `Annotations` can store the sequence and column mapping after operations like `gapstrip!` or `adjustreference!`.
  * Read and write `FASTA`, `Stockholm` or `Raw` formats.
  * Functions for shuffling the MSA: `shuffle_...`
* **Clustering** defines a `Clusters` type with the clusters information of the sequences in an MSA. `hobohmI` implements Hobohm I algorithm [4] and gives sequence weights according to the number of sequences in the clusters.
* **PDB** defines functions for parsing and working with `PDBFile` and `PDBML` formats:
  * Defines the types: `PDBResidue`, `PDBResidueIdentifier`, `PDBAtom`, `Coordinates`
  * Functions for estimation of `distance` , `contact` between `PDBResidue`s and type of interactions as: `vanderwaals`, `ionic`, etc.
  * Macros and functions for getting residues or atoms from a list/vector of residues, i.e.: `@residues`
* **SIFTS** has functions for downloading and parsing PDB SIFTS XML files:
  * `DataBase` and subtypes `db...` (i.e. `dbUniProt`) for a residue-level mapping between databases.
  * `siftsmapping` function to allow an easy-to-use `Dict` mapping.
  * `SIFTSResidue` is a immutable composite types with the SIFTS mapping for the residue.
* **Information** functions and types for measuring information content:
  * BLOSUM62 probabilities: `BLOSUM62_Pi` and `BLOSUM62_Pij`
  * `ResidueContingencyTables` for counting or store the probabilities of `Residues` (using the 20 residues BLOSUM62 alphabet) it also allows gaps count.
  * `InformationMeasure`s:
    * `Entropy`
    * `MutualInformation`
    * `MutualInformationOverEntropy`
  * Corrections to the MI values for co-evolution estimation:
    * `APC!` for MIp [5]
    * `buslje09` for the Z score and MIp from Buslje et. al. 2009 [1]
* **Pfam** has methods for working with *Pfam* alignments.
* **Utils** for common utils functions and types in MIToS.


References
==========

1. Buslje, C. M., Santos, J., Delfino, J. M., & Nielsen, M. **Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information.** *Bioinformatics 2009*, 25(9), 1125-1131.
2. Buslje, C. M., Teppa, E., Di Doménico, T., Delfino, J. M., & Nielsen, M. **Networks of high mutual information define the structural proximity of catalytic sites: implications for catalytic residue identification.** *PLoS Comput Biol 2010*, 6(11), e1000978-e1000978.
3. Altschul, S. F., Madden, T. L., Schäffer, A. A., Zhang, J., Zhang, Z., Miller, W., & Lipman, D. J. **Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.** *Nucleic acids research 1997*, 25(17), 3389-3402.
4. Hobohm, U., Scharf, M., Schneider, R., & Sander, C. **Selection of representative protein data sets.** *Protein Science 1992*, 1(3), 409-417.
5. Dunn, Stanley D., Lindi M. Wahl, and Gregory B. Gloor. **Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction.** *Bioinformatics 2008*, 24(3), 333-340.


Structural Bioinformatics Unit
==============================

[![FIL](http://mistic.leloir.org.ar/imgs/logo_horizontal.png)](http://www.leloir.org.ar/)
