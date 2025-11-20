# Model transformations between pangenome graphs
Whole Genome Alignment (WGA) graphs are bidirected graphs with nodes labeled with multiple sequence alignments of homologous genome fragments.
Similarly to WGAs, variation graphs (VG) are founded on bidirected graphs, but their nodes are labeled with single DNA sequences. Therefore, VGs do not represent full alignments, they only indicate related identical nucleotides.
## Transformation of WGAs into VGs
wga2vg is an algorithm that transforms a WGA representation of a sequence set into a compatible and compact VG representation. wga2vg requires as input a whole genome alignment graph in MAF format (Multiple Alignment File) and outputs a variation graph in GFA format (Graphical Fragment Assembly). WGA blocks consisting of a single sequence and built from uncovered fragments of the input sequences should be added to the input MAF file.

### Usage
1. Change directory:
```
cd wga2vg
```
2. Add missing blocks (i.e. blocks consisting of a single sequence) to your MAF file.
```
python addMissingBlocks.py <input_maf> <input_fasta> <output_maf>
```
4. Run wga2vg on MAF file that covers the entire input sequences (i.e. the output of the previous step):
```
python wga2vg.py <input_maf> <output_gfa>
```

## Transformation of WGAs into VGs
vg2wga and block-detector are two algorithms transforming VGs into WGAs. vg2wga requires as input a variation graph in GFA format (Graphical Fragment Assembly) and outputs a whole genome alignment graph in MAF format (Multiple Alignment File).
### vg2wga
In WGA graphs constructed using vg2wga, blocks correspond to the nodes of the input VGs.
#### Installation
To install vg2wga, download the repo, then run the following commands:
```
cd vg2wga
g++ -O3 vg2wga.cpp -o vg2wga
```
#### Usage
To run vg2wga, use the following command format:
```
./vg2wga <input_gfa> <output_maf>
```

### block-detector
block-detector is based on the algorithm of SibeliaZ-LCB, but operates on a variation graph instead of a compacted de Bruijn graph. The MSA in each block is performed using a modification of the library [poapy](https://github.com/ljdursi/poapy). block-detector requires as input a variation graph in GFA format (Graphical Fragment Assembly) and outputs a whole genome alignment graph in MAF format (Multiple Alignment File).
#### Usage
To run block-detector run the following-commands:
```
cd block_detector
python main.py -i <input_gfa> -o <output_maf>
```
The following not required flags can be used.
```
    -a - Abundance pruning parameter, default to 150. Vertices occurring more than a times are not considered as collinear walk seeds.
    -b - Maximal bubble size (for each walk), default to 200 residues.
    -m - Minimal collinear walk length, default to 50.
    --match - Match score in alignment (used if --align is True). Default to 5.
    --mismatch - Mismatch penalty in alignment (used if --align is True). Default to -4.
    --gap - Gap penalty in alignment (used if --align is True). Default to -8.
```
