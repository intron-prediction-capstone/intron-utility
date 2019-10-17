#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include "lib/gtf-cpp/gtf.h"
#include "lib/pdqsort/pdqsort.h"
#include "fasta.h"

const double EPSILON = 0.00000001111111111111;
#define DBL_EQ(a,b) (((a) <= (b) + EPSILON) && ((a) >= (b) - EPSILON))

// background frequency
const double BACKGROUND = 0.25;
const double ADJUSTMENT = 0.00001;

using Matrix = std::map<char, std::vector<double>>;

static Matrix pwm{
    {'A', {}},
    {'C', {}},
    {'T', {}},
    {'G', {}},
};

static void print_matrix(const Matrix& m) {
    for (auto& [key, val] : m) {
        std::cout << key;
        for (auto d : val) {
            std::cout << '\t' << d;
        }
        std::cout << '\n';
    }
}

static inline double score(double f) {
    if (DBL_EQ(f, 0.0)) return 0.0;
    double tmp = f / BACKGROUND;
    return log2(tmp);
}

static inline double normalize(double d, double min, double max) {
    if (DBL_EQ(d, 0.0)) {
        d = 50.0;
    } else if (d > 0.0) {
        // greater than or equal to 0 gets rescaled
        // to fit within 50-100
        // divide the value by the max and then multiply by 50, then
        // add 50 so that max is 100
        d = (((d * 100.0) / (max * 100.0)) * 50.0) + 50.0;
    } else {
        // negative values go from 0-50
        // take the value, divide by the min, then multiply by 50
        // subtract from 50 to get the opposite (we want the biggest value
        // to be at 0, because this is negative numbers)
        d = 50.0 - (((d * 100.0) / (min * 100.0)) * 50.0);
    }

    return d;
}

// TODO modify this
static Matrix normalize(const Matrix& m) {
    // find max and min values
    double min = 0.0,
           max = 0.0;
    for (auto& [key, val] : m) {
        for (double d : val) {
            if (d < min) min = d;
            if (d > max) max = d;
        }
    }
    
    // normalize
    Matrix ret = m;
    for (auto& [key, val] : ret) {
        for (double& d : val) {
            d = normalize(d, min, max);
        }
    }

    return ret;
}

int main(int argc, char** argv) {
    std::string pwmfilename, gtffilename, fastafilename;
    
    if (argc > 1) {
        pwmfilename = argv[1];
    }
    if (argc > 2) {
        gtffilename = argv[2];
    }
    if (argc > 3) {
        fastafilename = argv[3];
    }

    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s [pwm] [gtf] [fasta]\n", argv[0]);
        return 1;
    }

    std::ifstream pwmfile(pwmfilename);
    if (!pwmfile) {
        std::cerr << "Error opening file: " << pwmfilename << '\n';
        return 1;
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
                pwm[nt].push_back(tmp);
            }
        }
    }
    
    pwmfile.close();
    
    // open the GTF file
    GTFFile gtf;
    gtf.setfilename(gtffilename);
    try {
        gtf.load();
    } catch (const GTFError& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    std::vector<GTFSequence> exons = gtf.filter([](auto seq) -> bool {
        return seq.feature == "exon";
    });
    // sort exons by start
    pdqsort(exons.begin(), exons.end(), [](GTFSequence& a, GTFSequence& b) {
        return a.start < b.start;
    });

    std::cout << "Out of "
        << gtf.count()
        << " sequences, found "
        << exons.size()
        << " exons.\n";
    
    FASTAFile fasta;
    if (!fasta.open(fastafilename)) {
        std::cerr << "Error opening fasta file: " << fastafilename << '\n';
        return 1;
    }

    std::vector<std::string> introns;
    std::string tmpstr;
    for (std::size_t i = 0; i < exons.size() - 1; i++) {
        if (exons[i].end > exons[i+1].start) continue;
        try {
            tmpstr = fasta.get_sequence(exons[i].end - 2, exons[i].end + 10);
            tmpstr += " ... ";
            tmpstr += fasta.get_sequence(exons[i+1].start - 14, exons[i+1].start + 2);
            introns.push_back(tmpstr);
        } catch(const std::runtime_error& e) {
            std::cout << "Error: " << e.what() << '\n';
            std::cout << "i = " << i << '\n';
            std::cout << "start = " << exons[i].end - 2 << '\n'
                << "end = " << exons[i+1].start + 2 << '\n';
            return 1;
        }
    }
    introns.shrink_to_fit();

    for (std::size_t i = 0; i < 10; i++) {
        std::cout << introns[i] << '\n';
    }

    std::cout << "Loaded " << introns.size() << " introns.\n";

    return 0;
}
