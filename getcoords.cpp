#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <map>
#include <limits>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <utility>
#include <cmath>
#include <algorithm>
#include "lib/gtf-cpp/gtf.h"
#include "lib/pdqsort/pdqsort.h"
#include "lib/fasta-cpp/fasta.h"
#include "introns.h"

static int get_introns(const std::string& gtffile,
        std::vector<introns::Intron>& outputintrons) {

    // open the GTF file
    GTFFile gtf;
    gtf.setfilename(gtffile);
    try {
        gtf.load_filter([](const auto& seq) {
            return seq.feature == "exon";
        });
    } catch (const GTFError& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    std::vector<GTFSequence>& tmpexons = gtf.getall();

    // sort exons by start
    pdqsort(tmpexons.begin(), tmpexons.end(), [](GTFSequence& a, GTFSequence& b) {
        return a.start < b.start;
    });

    // split the exons into a map:
    /*
     * {
     *   gene_id: {
     *     transcript_id: {
     *       { exons, introns }
     *     }
     *   }
     * }
     */
    std::map<std::string, // gene_id
        std::map<std::string, // transcript_id
        std::pair<std::vector<GTFSequence>, // exons
        std::vector<introns::Intron>>>> // introns
            exons_by_transcript_by_gene;
    for (auto& exon : tmpexons) {
        exons_by_transcript_by_gene
            [exon.attributes["gene_id"]]
            [exon.attributes["transcript_id"]].first.push_back(exon);
    }

    // clean up memory as needed
    for (auto& [gene_id, tmap] : exons_by_transcript_by_gene) {
        for (auto& [tid, emap] : tmap) {
            emap.first.shrink_to_fit();
        }
    }

    std::cout << "-> Loaded " << tmpexons.size() << " exons\n";
    
    outputintrons.clear();
    std::string tmpstr;
    int tcount = 0;
    for (auto& [gene_id, exons_by_transcript] : exons_by_transcript_by_gene) {
        for (auto& [transcript_id, pair] : exons_by_transcript) {
            //if (tcount++ == 10) goto done;
            //std::cerr << "gene_id: " << gene_id << "; transcript " << transcript_id << '\n';
            auto& [exons, _introns] = pair;
            for (std::size_t i = 0; i < exons.size() - 1; i++) {
                // 5' is the first 14 chars, 3' is the last 18
                // the full sequence is everything minus the first 3 and last 4
                _introns.push_back({
                        exons[i].seqname,//sequence
                        "", // 3'
                        "",//tmpstr.substr(tmpstr.size() - 18), // 3'
                        "", // B'
                        "", // whole intron
                        exons[i].end-2, exons[i+1].start+3, // start and end
                        transcript_id,
                        gene_id,
                        { transcript_id }, // set the list of tids
                        .strand = exons[i].strand, // negative?
                        });
            }
        }
    }

done:

#if 0
    // keep a map of gene id to transcript id to duplicate counts
    std::map<std::string, std::map<std::string, unsigned long long>> dupecounts;

    for (auto& [gene_id, exons_by_transcript] : exons_by_transcript_by_gene) {
        for (auto& [transcript_id, pair] : exons_by_transcript) {
            auto& [exons, _introns] = pair;
            for (auto& intron : _introns) {
                for (auto& [ts_id, pair2] : exons_by_transcript) {
                    auto& [exons2, _introns2] = pair2;
                    for (auto& other : _introns2) {
                        if (&other == &intron) continue;
                        if (!other.keep_in_output) continue;
                        if (std::abs((long long)intron.start - (long long)other.start) <= 20L && std::abs((long long)intron.end - (long long)other.end) <= 20L) {
                            intron.keep_in_output = true;
                            dupecounts[gene_id][ts_id]++;
                            if (std::find(other.all_transcripts.begin(),
                                        other.all_transcripts.end(),
                                        transcript_id)
                                    == other.all_transcripts.end()) {
                                other.all_transcripts.push_back(transcript_id);
                            }
                        } else {
                            // if they overlap still, we need to check if one is
                            // really big and throw it out if it's huge
                            if (intron.start < other.end && intron.end > other.start) {
                                if (intron.full_sequence.size() > other.full_sequence.size()) {
                                    intron.keep_in_output = false;
                                } else if (other.full_sequence.size() > intron.full_sequence.size()) {
                                    other.keep_in_output = false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif

    for (auto& [gene_id, exons_by_transcript] : exons_by_transcript_by_gene) {
        for (auto& [transcript_id, pair] : exons_by_transcript) {
            auto& [exons, _introns] = pair;
            for (auto& intron : _introns) {
                if (intron.keep_in_output) {
                    outputintrons.push_back(intron);
                }
            }
        }
    }
            
    outputintrons.shrink_to_fit();

    std::cout << "-> Loaded " << outputintrons.size() << " introns\n";
    
#if 0
    std::time_t t = std::time(nullptr);
    char mbstr[100];
    std::strftime(mbstr, 100, "duplicatestats-%Y%m%d%H%M.csv", std::localtime(&t));

    std::ofstream outfile(mbstr);
    if (!outfile) {
        std::cerr << "Error opening " << mbstr << " for writing!\n"
            << "Continuing without duplicate stats output.\n";
    } else {
        outfile << "\"gene_id\",\"transcript_id\",\"duplicates\"\r\n";
        for (auto& [gid, ts] : dupecounts) {
            for (auto& [tid, dupes] : ts) {
                outfile << '"' << gid << "\",\"" << tid << "\"," << dupes << "\r\n";
            }
        }
        outfile.close();
        std::cout << "-> Wrote duplicate stats to " << mbstr << '\n';
    }
#endif

    return 0;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr <<
R"usage(
usage: getcoords <gtf> <out+> <out->
     gtf:  the GTF file to read from
     out+: the file to write positive strand introns to
     out-: the file to write negative strand introns to
)usage";
        return 1;
    }
    std::string fname = argv[1];
    std::string posfile = argv[2];
    std::string negfile = argv[3];

    std::vector<introns::Intron> introns;
    get_introns(fname, introns);

    std::ofstream outpos(posfile);
    if (!outpos) {
        std::cerr << "Could not open positive file for writing!\n";
        return 1;
    }

    std::ofstream outneg(negfile);
    if (!outneg) {
        std::cerr << "Could not open negative file for writing!\n";
        return 1;
    }

    for (auto& i : introns) {
        (i.strand == '-' ? outneg : outpos) << i.gene_id << ',' << i.transcript_id << ',' << i.sequence << ':' << i.start << '-' << i.end << '\n';
    }

    return 0;
}
