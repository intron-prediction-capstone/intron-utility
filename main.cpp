#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <map>
#include <limits>
#include <cstdio>
#include <cstring>
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

static void usage(char* executable) {
    std::fprintf(stderr,
            "Usage:\n"
            "\t1: %s --parse [gtf] [fasta] [output]\n"
            "\t2: %s [3' pwm] [5' pwm] [B' pwm] [gtf] [fasta]\n"
            "\n\t\t1: Parse a FASTA file with a GTF file and print out any introns"
            "\n\t\t   found into the [output] file."
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
    if (argc > 4) {
        if (parseronly) {
            outputfilename = argv[4];
        } else {
            gtffilename = argv[4];
        }
    }
    if (!parseronly) {
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
        std::ofstream outfile(outputfilename);
        if (!outfile) {
            std::cerr << "Error opening output file for writing: " << outputfilename << '\n';
            return 1;
        }

        int ret = get_introns(gtffilename, fastafilename, introns);
        if (ret == 0) {
            std::cout << "-> Writing introns to file...";
            std::size_t intronswritten = 0;
            const char* logfmt = "\r-> Writing introns to file: [%lu/%lu]";
            for (auto& intron : introns) {
                outfile << intron.full_sequence << '\n';
                intronswritten++;
                std::printf(logfmt, intronswritten, introns.size());
            }
            std::cout << "\n-> Wrote introns to file\n";
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
        gtf.load();
    } catch (const GTFError& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    std::vector<GTFSequence> tmpexons = gtf.filter([](const auto& seq) {
        return seq.feature == "exon";
    });

    // sort exons by start
    pdqsort(tmpexons.begin(), tmpexons.end(), [](GTFSequence& a, GTFSequence& b) {
        return a.start < b.start;
    });

    // split the exons into a map of gene id -> exons
    std::map<std::string, std::vector<GTFSequence>> exons_by_gene;
    for (auto& exon : tmpexons) {
        try {
            exons_by_gene.at(exon.attributes["gene_id"]).push_back(exon);
        } catch(const std::out_of_range& oor) {
            exons_by_gene.insert({exon.attributes["gene_id"], std::vector<GTFSequence>(1, exon)});
        }
    }

    // clean up memory as needed
    for (auto& [key, val] : exons_by_gene) {
        val.shrink_to_fit();
    }

    std::cout << "-> Loaded " << tmpexons.size() << " exons\n";
    
    FASTAFile fasta;
    if (!fasta.open(fastafile)) {
        std::cerr << "Error opening fasta file: " << fastafile << '\n';
        return 1;
    }

    outputintrons.clear();
    std::string tmpstr;
    for (auto& [gene_id, exons] : exons_by_gene) {
        for (std::size_t i = 0; i < exons.size() - 1; i++) {
            // cut out overlapped exons
            if (exons[i].end > exons[i+1].start) continue;
            try {
                // 3nt from upstream
                // 4nt from downstream
                tmpstr = fasta.get_sequence(exons[i].end - 2, exons[i+1].start + 3, true);

                // 5' is the first 14 chars, 3' is the last 18
                // the full sequence is everything minus the first 3 and last 4
                outputintrons.push_back({
                        tmpstr.substr(0, 14), // 5'
                        tmpstr.substr(tmpstr.size() - 18), // 3'
                        tmpstr.substr(3, tmpstr.size() - 7), // whole intron
                        exons[i].end+1, exons[i+1].start-1, // start and end
                        exons[i].attributes["transcript_id"], // gene id
                        gene_id,
                    });
            } catch(const std::runtime_error& e) {
                std::cerr << "Error: " << e.what() << '\n';
                return 1;
            }
        }
    }
    outputintrons.shrink_to_fit();

    std::cout << "-> Loaded " << outputintrons.size() << " introns\n";

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
