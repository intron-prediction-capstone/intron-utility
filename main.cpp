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


// 3' PWM
static introns::Matrix pwm_lod_3p;

// 5' PWM
static introns::Matrix pwm_lod_5p;

// B' PWM
static introns::Matrix pwm_lod_Bp;

// parse a pwm file and return the pwm
// this also applies the `score' function to the values to make a lod
static introns::Matrix parse_pwm(const std::string&);

// parser-only operation
static int get_introns(const std::string& gtffile,
        const std::string& fastafile,
        std::vector<introns::Intron>& outputintrons);

// scores and normalizes introns
static void score_and_normalize_introns(std::vector<introns::Intron>& introns);

static void output_introns_fasta(std::vector<introns::Intron>& introns);
static void output_introns_gtf(std::vector<introns::Intron>& introns);

static void usage(char* executable) {
    std::fprintf(stderr,
            "Usage:\n"
            "\t1: %s --parse [gtf] [fasta]\n"
            "\t2: %s [3' pwm] [5' pwm] [B' pwm] [gtf] [fasta]\n"
            "\n\t\t1: Parse a FASTA file with a GTF file.\n"
            "\n\t\t2: Take a 3' PWM, 5' PWM, and B' PWM and predict U12 introns from"
            "\n\t\t   the given GTF and FASTA files.\n",
            executable, executable);
}

int main(int argc, char** argv) {
    std::string pwmfilename_3p,
        pwmfilename_5p,
        pwmfilename_Bp,
        gtffilename,
        fastafilename,
        outputfilename;

    bool parseronly = false;

    if (argc == 1) {
        usage(argv[0]);
        return 1;
    }

    if (argc > 1) {
        parseronly = (0 == std::strcmp(argv[1], "--parse"));
        if (!parseronly) {
            pwmfilename_3p = argv[1];
        }
    }
    if (argc > 2) {
        if (parseronly) {
            gtffilename = argv[2];
        } else {
            pwmfilename_5p = argv[2];
        }
    }
    if (argc > 3) {
        if (parseronly) {
            fastafilename = argv[3];
        } else {
            pwmfilename_Bp = argv[3];
        }
    }
    if (!parseronly) {
        if (argc > 4) {
            gtffilename = argv[4];
        }
        if (argc > 5) {
            fastafilename = argv[5];
        }
    }

    if ((parseronly && argc < 5) || (!parseronly && argc < 6)) {
        usage(argv[0]);
        return 1;
    }

    // do parser-only operation
    if (parseronly) {
        std::vector<introns::Intron> introns;

        int ret = get_introns(gtffilename, fastafilename, introns);
        if (ret == 0) {
            output_introns_fasta(introns);
            output_introns_gtf(introns);
        }

        return 0;
    }

    pwm_lod_3p = parse_pwm(pwmfilename_3p);
    pwm_lod_5p = parse_pwm(pwmfilename_5p);
    pwm_lod_Bp = parse_pwm(pwmfilename_Bp);
    std::cout << "-> Loaded PWMs\n";
    
    std::vector<introns::Intron> introns;
    // get the introns
    int ret = get_introns(gtffilename, fastafilename, introns);
    if (ret != 0) return ret;

    score_and_normalize_introns(introns);

    std::cout << "-> Scored introns\n";

    output_introns_fasta(introns);
    output_introns_gtf(introns);
    
    return 0;
}

// get introns
static int get_introns(const std::string& gtffile,
        const std::string& fastafile,
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
    
    FASTAFile fasta;
    if (!fasta.open(fastafile)) {
        std::cerr << "Error opening fasta file: " << fastafile << '\n';
        return 1;
    }

    outputintrons.clear();
    std::string tmpstr;
    int tcount = 0;
    for (auto& [gene_id, exons_by_transcript] : exons_by_transcript_by_gene) {
        for (auto& [transcript_id, pair] : exons_by_transcript) {
            if (tcount++ == 10) goto done;
            std::cerr << "gene_id: " << gene_id << "; transcript " << transcript_id << '\n';
            auto& [exons, _introns] = pair;
            for (std::size_t i = 0; i < exons.size() - 1; i++) {
                // cut out overlapped exons
                if (exons[i].end > exons[i+1].start) continue;
                try {
                    // 3nt from upstream
                    // 4nt from downstream
                    tmpstr = fasta.get_sequence(exons[i].end - 2, exons[i+1].start + 3, true);

                    // 5' is the first 14 chars, 3' is the last 18
                    // the full sequence is everything minus the first 3 and last 4
                    _introns.push_back({
                            tmpstr.substr(0, 14), // 5'
                            tmpstr.substr(tmpstr.size() - 18), // 3'
                            tmpstr.substr(3, tmpstr.size() - 7), // whole intron
                            exons[i].end+1, exons[i+1].start-1, // start and end
                            transcript_id,
                            gene_id,
                            { transcript_id }, // set the list of tids
                            .strand = exons[i].strand, // negative?
                            });

                    // if this is negative, reverse + complement the strands
                    if (_introns.back().strand == '-') {
                        _introns.back() = !(_introns.back());
                    }
                } catch(const std::runtime_error& e) {
                    std::cerr << "Error: " << e.what() << '\n';
                    return 1;
                }
            }
        }
    }

done:

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

    return 0;
}

// scores and normalizes introns
static void score_and_normalize_introns(std::vector<introns::Intron>& introns) {
    for (auto& intron : introns) {
        double sum = 0.0;
        for (std::size_t i = 0; i < intron.five_prime.size(); i++) {
            sum += pwm_lod_5p[intron.five_prime[i]][i];
        }
        for (std::size_t i = 0; i < intron.three_prime.size(); i++) {
            sum += pwm_lod_3p[intron.three_prime[i]][i];
        }
        intron.score = sum;
    }

    double min = std::numeric_limits<double>::max(),
           max = std::numeric_limits<double>::lowest();
    for (auto& i : introns) {
        if (i.score > max) max = i.score;
        if (i.score < min) min = i.score;
    }

    for (auto& i : introns) {
        i.score_normalized = introns::normalize(i.score, min, max);
    }
}

// parse a pwm from a file
// this also applies the `score' function to the values to make a lod
static introns::Matrix parse_pwm(const std::string& pwmfilename) {
    introns::Matrix _pwm{
        {'G', {}},
        {'T', {}},
        {'A', {}},
        {'C', {}},
    };

    std::ifstream pwmfile(pwmfilename);
    if (!pwmfile) {
        throw std::runtime_error("Error opening file: " + pwmfilename);
    }

    char nt;
    std::stringstream ss;
    std::string line;
    double tmp;
    int count;
    while (std::getline(pwmfile, line)) {
    	count = 0;
        ss = std::stringstream(line);
        while (ss) {
            ss >> tmp;
            if (ss) {
            	count++;
            	if(count == 1) {
            		nt = 'A';
            	} else if (count == 2) {
            		nt = 'C';
            	} else if (count == 3) {
            		nt = 'G';
            	} else if (count == 4) {
            		nt = 'T';
            	}
                _pwm[nt].push_back(introns::log_odds(tmp));
            }
        }
    }
    
    pwmfile.close();

    return _pwm;
}

static void output_introns_fasta(std::vector<introns::Intron>& introns) {
    std::time_t t = std::time(nullptr);
    char mbstr[100];
    std::strftime(mbstr, 100, "introns-%Y%m%d%H%M.fa", std::localtime(&t));

    std::ofstream outfile(mbstr);
    if (!outfile) {
        throw std::runtime_error("Error opening file for writing: " + std::string(mbstr));
        return;
    }

    std::cout << "-> Writing sequences to " << mbstr << '\n';

    for (auto& intron : introns) {
        outfile << '>' << intron.gene_id
            << ',' << intron.transcript_id
            << ',' << intron.start
            << '-' << intron.end << "\r\n";
        // print the intron 60 characters at a time
        std::size_t start = 0;
        std::string tmp = "";
        while (start < intron.full_sequence.size()) {
            tmp = intron.full_sequence.substr(start, 60);
            outfile << tmp << "\r\n";
            start += 60;
        }
    }

    outfile.close();
}

static void output_introns_gtf(std::vector<introns::Intron>& introns) {
    std::time_t t = std::time(nullptr);
    char mbstr[100];
    std::strftime(mbstr, 100, "introns-%Y%m%d%H%M.gtf", std::localtime(&t));

    GTFFile outfile(mbstr);
    if (!outfile.open_for_writing(false)) {
        throw std::runtime_error("Error opening file for writing: " + std::string(mbstr));
        return;
    }

    std::cout << "-> Writing sequence annotations to " << mbstr << '\n';

    std::size_t id = 0;
    GTFSequence tmpseq;
    for (auto& intron : introns) {
        tmpseq = {
            std::to_string(id++),
            "intron-util",
            "intron_CNS",
            intron.start,
            intron.end,
            (DBL_EQ(intron.score_normalized, 0.0) ? '.' : intron.score_normalized),
            intron.strand,
            -1,
            {
                { "gene_id", intron.gene_id },
                { "transcript_id", intron.transcript_id },
            },
        };
        outfile << tmpseq;
    }

    outfile.close();
}
