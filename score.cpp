#include "introns.h"
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cctype>
#include <stdexcept>
#include <utility>
#include "lib/gtf-cpp/gtf.h"

// lines in the list files which tell us things about each intron we parse
struct listline {
    std::string geneid, transcriptid, sequence;
    std::size_t start, end;
};

// chunk of data from the FASTA file
struct parsedchunk {
    std::string title, data;
};

// 3' PWM
static introns::Matrix pwm_lod_3p;

// 5' PWM
static introns::Matrix pwm_lod_5p;

// B' PWM
static introns::Matrix pwm_lod_Bp;

static double score_bp(const std::string& bp) {
    double s = 0.0;
    for (std::size_t i = 0; i < 12; i++) {
        if (bp[i] == 'A' || bp[i] == 'C'
                || bp[i] == 'T' || bp[i] == 'G') {
            s += pwm_lod_Bp[bp[i]][i];
        }
    }
    return s;
}

// returns a pair of {index, score} for the branch point
static std::pair<std::size_t, double> find_bp(const std::string& intron) {
    std::size_t maxoff = 0;
    double maxscore = 0.0;
    double s;
    for (std::size_t i = 40; i >= 1; i--) {
        s = score_bp(intron.substr(intron.size()-14-i, 12));
        if (s > maxscore) {
            maxscore = s;
            maxoff = intron.size()-14-i;
        }
    }
    return std::make_pair(maxoff, maxscore);
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

// scores introns
static void score_and_normalize_introns(std::vector<introns::Intron>& introns) {
    for (auto& intron : introns) {
        intron.valid = true;
        if (intron.full_sequence.length() <= 54) {
            intron.valid = false;
            continue;
        }
        double sum = 0.0;
        for (std::size_t i = 0; i < intron.five_prime.size(); i++) {
            double s = pwm_lod_5p[intron.five_prime[i]][i];
            intron.fivescore = s;
            sum += s;
        }
        for (std::size_t i = 0; i < intron.three_prime.size(); i++) {
            double s = pwm_lod_3p[intron.three_prime[i]][i];
            intron.threescore = s;
            sum += s;
        }
        intron.score = sum;

        auto [bpi, bps] = find_bp(intron.full_sequence);
        intron.bpscore = bps;
        intron.score += bps;
        intron.branch_point = intron.full_sequence.substr(bpi, 12);
        intron.bp_index = bpi;
    }

    double min = std::numeric_limits<double>::max(),
           max = std::numeric_limits<double>::lowest();
    for (auto& i : introns) {
        if (!i.valid) continue;
        if (i.score > max) max = i.score;
        if (i.score < min) min = i.score;
    }

    for (auto& i : introns) {
        if (!i.valid) {
            i.score_normalized = -1.0;
        } else {
            i.score_normalized = introns::normalize(i.score, min, max);
        }
    }
}

// load and score introns; can be run in a separate thread
// If an error happens, it will print and then return an empty vector
std::vector<introns::Intron> load_and_parse(std::ifstream& inlist,
        std::ifstream& inparsed, bool isnegative) {
    // for each line in the list, we should have one chunk in the parsed
    // file (FASTA)
    std::vector<listline> list;
    listline tmpline;
    std::string tmpstr;
    while (std::getline(inlist, tmpline.geneid, ',')) {
        std::getline(inlist, tmpline.transcriptid, ',');
        std::getline(inlist, tmpline.sequence, ':');
        std::getline(inlist, tmpstr, '-');
        tmpline.start = (std::size_t)std::atoll(tmpstr.c_str());
        std::getline(inlist, tmpstr);
        tmpline.end = (std::size_t)std::atoll(tmpstr.c_str());
        list.push_back(tmpline);
    }
    list.shrink_to_fit();

    std::printf("Loaded %zu %s strands' info\n", list.size(),
            isnegative? "negative":"positive");

    std::vector<parsedchunk> chunks;
    //chunks.reserve(list.size()); // don't reallocate this vector to be huge
    // we have read every item in the list, now we can read everything in
    // the FASTA
    char c;
    bool hasfirst = false;
    parsedchunk tmpchunk;
    while ((c = inparsed.get()) != std::ifstream::traits_type::eof()) {
        if (c == '>') {
            if (!hasfirst) {
                hasfirst = true;
            } else {
                chunks.push_back(tmpchunk);
            }
            tmpchunk.data.clear();
            std::getline(inparsed, tmpchunk.title);
        } else {
            if (std::isalpha(c)) {
                tmpchunk.data += c;
            }
        }
    }
    chunks.push_back(tmpchunk);
    chunks.shrink_to_fit();

    std::printf("Loaded %zu %s strand chunks\n", chunks.size(),
            isnegative? "negative":"positive");

    std::vector<introns::Intron> intr;
    //intr.reserve(list.size());
    introns::Intron tmpintron;
    for (std::size_t i = 0; i < list.size(); i++) {
        try {
            tmpintron = {
                list[i].sequence,
                chunks[i].data.substr(0, 14),
                chunks[i].data.substr(chunks[i].data.size() - 18),
                "",
                chunks[i].data.substr(3, chunks[i].data.size() - 7),
                list[i].start + 3,
                list[i].end - 4,
                list[i].transcriptid,
                list[i].geneid,
                { list[i].transcriptid },
            };
            tmpintron.strand = isnegative;
            intr.push_back(tmpintron);
        } catch (std::out_of_range& e) {
            continue;
        }
    }
    intr.shrink_to_fit();

    return intr;
}

static void output_introns_fasta(std::vector<introns::Intron>& introns);
static void output_introns_gtf(std::vector<introns::Intron>& introns);

int main(int argc, char** argv) {
    if (argc < 8) {
        std::cerr << "usage: score <3'pwm> <5'pwm> <B'pwm> <poslist> <neglist> <posparsed> <negparsed>\n";
        return 1;
    }
    std::string pwmfilename_3p,
        pwmfilename_5p,
        pwmfilename_Bp;
    
    pwmfilename_3p = argv[1];
    pwmfilename_5p = argv[2];
    pwmfilename_Bp = argv[3];

    pwm_lod_3p = parse_pwm(pwmfilename_3p);
    pwm_lod_5p = parse_pwm(pwmfilename_5p);
    pwm_lod_Bp = parse_pwm(pwmfilename_Bp);

    std::string poslist = argv[4];
    std::string neglist = argv[5];
    std::string posparsed = argv[6];
    std::string negparsed = argv[7];

    std::vector<introns::Intron> allintrons;

    std::ifstream inposlist(poslist);
    if (!inposlist) {
        std::cerr << "Could not open file " << poslist << '\n';
        return 1;
    }
    std::ifstream inneglist(neglist);
    if (!inneglist) {
        std::cerr << "Could not open file " << neglist << '\n';
        return 1;
    }
    std::ifstream inposparsed(posparsed);
    if (!inposparsed) {
        std::cerr << "Could not open file " << posparsed << '\n';
        return 1;
    }
    std::ifstream innegparsed(negparsed);
    if (!innegparsed) {
        std::cerr << "Could not open file " << negparsed << '\n';
        return 1;
    }

    auto posintrons = load_and_parse(inposlist, inposparsed, false);
    auto negintrons = load_and_parse(inneglist, innegparsed, true);

    allintrons = posintrons;
    allintrons.insert(allintrons.end(), negintrons.begin(), negintrons.end());

    std::cout << "Got results; scoring...\n";

    score_and_normalize_introns(allintrons);

    std::cout << "Scoring complete; writing output...\n";

    // output data
    output_introns_fasta(allintrons);
    output_introns_gtf(allintrons);
    return 0;
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

    GTFSequence tmpseq;
    for (auto& intron : introns) {
        tmpseq = {
            intron.sequence,
            "intron-util",
            "intron_CNS",
            intron.start,
            intron.end,
            (DBL_EQ(intron.score_normalized, 0.0) ? '.' : intron.score_normalized),
            intron.strand ? '-' : '+',
            -1,
            {
                { "gene_id", intron.gene_id },
                { "transcript_id", intron.transcript_id },
                { "five_score", std::to_string(intron.fivescore) },
                { "three_score", std::to_string(intron.threescore) },
                { "bp_score", std::to_string(intron.bpscore) },
                { "bp", intron.branch_point },
            },
        };
        outfile << tmpseq;
    }

    outfile.close();
}
