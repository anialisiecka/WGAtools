# Model transformations

## Transformation of WGAs into VGs
wga2vg is an algorithm that transforms a WGA representation of a sequence set into a compatible and compact VG representation.

### Usage
1. Change directory:
```
cd wga2vg
```
2. Add missing blocks (i.e. blocks consisting of a single sequence) to your MAF file.
```
python addMissingBlocks.py <input_maf> <input_fasta> <output_maf>
```
4. Run wga2vg:
```
python wga2vg.py <input_maf> <output_gfa>
```

## Transformation of WGAs into VGs
vg2wga and block-detector are two algorithms transforming VGs into WGAs.
In WGA graphs constructed using vg2wga, blocks correspond to the nodes of the input VGs.
### Installation of vg2wga
To install vg2wga, download the repo, then run the following commands:
```
cd vg2wga
g++ -O3 vg2wga.cpp -o vg2wga
```
### Usage of vg2wga
To run vg2wga, use the following command format:
```
./vg2wga <input_gfa> <output_maf>
```
block-detector is based on the algorithm of SibeliaZ-LCB, but operates on a variation graph instead of a compacted de Bruijn graph. The MSA in each block is performed using a modification of the library poapy.

## Usage of block-detector
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

## License
This project is licensed under the MIT License - see the [LICENSE](./LICENSE) file for details.
