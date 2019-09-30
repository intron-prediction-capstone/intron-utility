#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#define EPSILON 0.00000001111111111111
#define DBL_EQ(a,b) (((a) <= (b) + EPSILON) && ((a) >= (b) - EPSILON))

// takes char 'n' and makes sure that it's a capital letter
// this means we can make 'a' -> 'A' and so on to keep things consistent
#define CAP_NT(n) ((n)>='a'&&(n)<='z'?(n)+('A'-'a'):(n))

// background frequency
#define BACKGROUND 0.25
#define ADJUSTMENT 0.00001

#define INPUT_FILE "input.tsv"

using Matrix = std::map<char, std::vector<double>>;

static Matrix frequencies;
static Matrix pwm;

static void print_matrix(const Matrix& m) {
    for (auto& it : m) {
        std::cout << it.first;
        for (auto d : it.second) {
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

static Matrix normalize(const Matrix& m) {
    // find max and min values
    double min = 0.0,
           max = 0.0;
    for (auto& pair : m) {
        for (double d : pair.second) {
            if (d < min) min = d;
            if (d > max) max = d;
        }
    }
    
    // normalize
    Matrix ret = m;
    for (auto& pair : ret) {
        for (double& d : pair.second) {
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
        }
    }

    return ret;
}

int main(void) {
    std::ifstream infile(INPUT_FILE);
    if (!infile) {
        std::cerr << "Error reading input file!\n";
        return 1;
    }

    char nt;
    std::vector<double> freq;
    std::stringstream ss;
    std::string line;
    double tmp;
    while (std::getline(infile, line)) {
        ss = std::stringstream(line);
        ss >> nt;
        while (ss) {
            ss >> tmp;
            if (ss) {
                freq.push_back(tmp);
            }
        }
        frequencies[CAP_NT(nt)] = freq;
        freq.clear();
    }
    
    infile.close();

    // confirm that all lines are the same length
    std::size_t freq_count = frequencies['A'].size();
    for (auto& pair : frequencies) {
        if (pair.second.size() != freq_count) {
            std::cerr << "Error: not all lines are the same length!\n";
            return 1;
        }
    }

    std::cout << "Input Matrix:\n";
    print_matrix(frequencies);

    std::vector<double> scores;
    for (auto& pair : frequencies) {
        pwm[pair.first] = std::vector<double>(freq_count, 0.0);
        for (std::size_t i = 0; i < freq_count; i++) {
            pwm[pair.first][i] = score(pair.second[i]);
        }
    }

    std::cout << "\nPWM Scores:\n" << std::fixed;
    std::cout.precision(3);
    print_matrix(pwm);

    std::cout << "\nNormalized PWM Scores:\n";
    std::cout.precision(1);
    print_matrix(normalize(pwm));

    return 0;
}
