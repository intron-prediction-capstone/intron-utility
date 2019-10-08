#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <iostream>

class FASTAFile {
    public:
        FASTAFile(): file("") {}
        FASTAFile(const std::string& filename): file(filename) {
            if (!open(filename)) {
                throw std::runtime_error("Error opening file: " + filename + "!");
            }
        }

        ~FASTAFile() { infile.close(); }

        // returns false if opening the file failed
        // use this only if you used the default constructor.
        bool open(const std::string& filename) {
            file = filename;
            infile = std::ifstream(file);
            // now we need to chop the first line out and figure out how long it is
            std::string tmp;
            std::getline(infile, tmp);
            header_len = tmp.size() + 1; // +1 for newline character
            // now get the length of each line for math later on
            std::getline(infile, tmp);
            line_nt = tmp.size(); // this is length in nt, not counting newline
            return (bool)infile;
        }

        // closes the file
        void close() {
            infile.close();
        }

        // gets a string of nucleotides from start to end, inclusive;
        // so specifying 1, 2 would get 2nt
        std::string get_sequence(std::size_t start, std::size_t end) {
            std::string ret;
            char tmp;
            infile.seekg(seq_start(start));
            std::size_t count = 0;
            while (count < (end - start) + 1) {
                tmp = infile.get();

                if (tmp == std::ifstream::traits_type::eof()) {
                    throw std::runtime_error("End coordinate out of bounds");
                }

                if (nullptr != std::strchr("actgACTGNn", tmp)) {
                    ret += tmp;
                    count++;
                }
            }
            return ret;
        }

    private:
        std::size_t header_len;
        std::size_t line_nt;
        std::string file;
        std::ifstream infile;

        // get the *actual* starting byte to seek to
        inline std::size_t seq_start(std::size_t start) {
            return (header_len + ((start / line_nt) * (line_nt + 1)) + (start % line_nt)) - 1;
        }
};

#endif /* FASTA_H */
