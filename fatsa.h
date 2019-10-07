#ifndef FATSA_H
#define FATSA_H

#include <fstream>
#include <string>
#include <stdexcept>

class FATSAFile {
    public:
        FATSAFile(): file("") {}
        FATSAFile(const std::string& filename): file(filename) {
            if (!open(filename)) {
                throw std::runtime_error("Error opening file: " + filename + "!");
            }
        }

        ~FATSAFile() { infile.close(); }

        bool open(const std::string& filename) {
            file = filename;
            infile = std::ifstream(file);
            // now we need to chop the first line out and figure out how long it is
            std::string tmp;
            infile.getline(tmp);
            header_len = tmp.size() + 1; // +1 for newline character
            // now get the length of each line for math later on
            infile.getline(tmp);
            line_nt = tmp.size(); // this is length in nt, not counting newline
            return (bool)infile;
        }

        void close() {
            infile.close();
        }

        // gets a string of nucleotides from start to end, inclusive;
        // so specifying 1, 2 would get 2nt
        std::string get_sequence(std::size_t start, std::size_t end) {
            std::string ret;
            infile.seekg(seq_start(start));
            // TODO
        }

    private:
        std::size_t header_len;
        std::size_t line_nt;
        std::string file;
        std::ifstream infile;

        // get the *actual* starting byte to seek to
        inline std::size_t seq_start(std::size_t start) {
            // TODO validate that this works as intended
            return header_len + ((start / line_nt) * (line_nt + 1)) - 1;
        }
};

#endif /* FATSA_H */
