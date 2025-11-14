#include <fstream>
#include <string>
#include <sstream>
#include <vector> 
#include <iostream>
#include <map> 
#include <utility>
#include <cmath>
#include <algorithm>
#include <set>

class Block {
        public:
                std::string genome;
                char strand;
                int start;
};

int main(int argc, char* argv[])
{
        std::string input_name = argv[1];
        std::string output_name = argv[2];
        std::ifstream gfa_file(input_name);
        std::string line;
        std::vector<std::string> vtxs;
        vtxs.push_back("");
        std::map<std::string, int> srcSizes;
        while (std::getline(gfa_file, line))
        {
                if (line[0] == 'S')
                {
                        std::stringstream ss(line);
                        std::string tag, vtx, sequence;
                        ss >> tag >> vtx >> sequence;
                        vtxs.push_back(sequence);
                } 
        }
        std::vector<std::vector<Block>> blocks(vtxs.size());
        std::ifstream gfa_file2(input_name);

        while (std::getline(gfa_file2, line))
        {
                if (line[0] == 'P')
                {
                        int srcSize = 0;
                        std::stringstream ss(line);
                        std::string tag, genome, path, cigar;
                        ss >> tag >> genome >> path >> cigar;
                        std::stringstream sp(path);
                        std::string oriented_vtx;

                        int pos = 0;
                        while (std::getline(sp, oriented_vtx, ','))
                        {
                                Block block;
                                char strand = oriented_vtx.back();
                                block.strand = strand;
                                block.genome = genome;
                                oriented_vtx.pop_back();
                                int vtx_id = std::stoi(oriented_vtx);
                                srcSize += vtxs[vtx_id].size();
                                int start;
                                if (strand == '-')
                                {
                                        start =  pos +  vtxs[vtx_id].size();
                                } else
                                {
                                        start = pos;
                                }
                                block.start = start;
                                pos += vtxs[vtx_id].size();
                                blocks[vtx_id].push_back(block);
                        }
                        srcSizes[genome] = srcSize;
                }
        }
        std::ofstream maf_file(output_name);
        maf_file << "##maf version=1 scoring=tba.v8\n\n";
        for (int i = 1; i < blocks.size(); i++)
        {
                if ((blocks[i].size()<=45)&&(blocks[i].size()>1)) {
                        maf_file << "a\n";
                        for (auto & block: blocks[i])
                        {
                                if (block.strand == '-')
                                {
                                        block.start = srcSizes[block.genome] - block.start;
                                }
                                maf_file << "s\t" << block.genome << "\t" << block.start << "\t" << vtxs[i].size() << "\t" << block.strand << "\t" << srcSizes[block.genome] << "\t" << vtxs[i] << "\n";
                        }
                        maf_file << "\n";
                }
        }
        maf_file.close();
        return 0;

}
