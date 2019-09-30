#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#define EPSILON 0.00000001111111111111
#define DBL_EQ(a,b) ((a <= b + EPSILON) && (a >= b - EPSILON))

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

int main(void) {
    std::ifstream infile(INPUT_FILE);
    if (!infile) {
        std::cout << "Error reading input file!\n";
        return 1;
    }

    char nt;
    std::vector<double> freq(5, 0.0);
    for (int i = 0; i < 4; i++) {
        infile >> nt;
        for (int j = 0; j < 5; j++) {
            infile >> freq[j];
        }
        frequencies[nt] = freq;
    }

    infile.close();

    std::cout << "Input Matrix:\n";
    print_matrix(frequencies);

    std::vector<double> scores;
    for (auto& pair : frequencies) {
        pwm[pair.first] = std::vector<double>(5, 0.0);
        for (int i = 0; i < 5; i++) {
            pwm[pair.first][i] = score(pair.second[i]);
        }
    }

    std::cout << "\nPWM Scores:\n" << std::fixed;
    std::cout.precision(3);
    print_matrix(pwm);

    return 0;
}
