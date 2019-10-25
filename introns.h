#ifndef INTRONS_H
#define INTRONS_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>

#define DBL_EQ(a,b) (((a) <= (b) + EPSILON) && ((a) >= (b) - EPSILON))
const double EPSILON = 0.00000001111111111111;

// background frequency
const double BACKGROUND = 0.25;
const double ADJUSTMENT = 0.00001;


namespace introns {
    using Matrix = std::map<char, std::vector<double>>;

    // Structure representing an intron
    struct Intron {
        // 5' sequence
        std::string five_prime;
        // 3' sequence
        std::string three_prime;
        // full intron *not counting upstream and downstream exon nt*
        std::string full_sequence;
        // start and end point (inclusive; no exonic nts)
        std::size_t start, end;
        // transcript and gene ID
        std::string transcript_id, gene_id;
        // raw score and normalized score
        double score, score_normalized;
        // after processing, this will be used to store all transcripts in which
        // this intron can be found
        std::vector<std::string> all_transcripts;
        // whether or not we want to keep this in output
        bool keep_in_output = true;
    };

    // prints a matrix
    inline void print_matrix(const introns::Matrix& m) {
        for (auto& [key, val] : m) {
            std::cout << key;
            for (auto d : val) {
                std::cout << '\t' << d;
            }
            std::cout << '\n';
        }
    }

    // calculates the log-odds score for a value
    inline double log_odds(double f) {
        if (DBL_EQ(f, 0.0)) return 0.0;
        return log2(f / BACKGROUND + ADJUSTMENT); // adjustment so that no log2(0)
    }

    // normalizes a value
    inline double normalize(double d, double min, double max) {
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
};

#endif /* INTRONS_H */
