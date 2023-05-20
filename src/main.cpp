#include <iostream>

#include "ReadMapping.h"
#include "logging.h"

using namespace std;

int main(int argc, char *argv[]) {
    // Parse command line arguments
    bool map_mode = false;
    string output_file = "output/mappings.txt";

    string reference_file, input_file;
    if (argc < 3 || argc > 5) {
        cout << "Usage: ./release reference_file input_file [output_file] [-m]" << endl;
        for (int i = 0; i < argc; i++) {
            cout << argv[i] << endl;
        }
        return 1;
    }

    reference_file = argv[1];
    input_file = argv[2];
    for (int i = 3; i < argc; i++) {
        if (string(argv[i]) == "-m") {
            map_mode = true;
        } else {
            output_file = argv[i];
        }
    }

    auto logger = new Logger(LogLevel::INFO_ONLY);

    logger->log_msg(LogLevel::INFO, "main", "Reading input files...");
    ReadMapping read_mapping(reference_file, input_file, logger);
    read_mapping.map_mode = map_mode;

//    vector<Mapping *> mappings;
//    read_mapping.map_read(new ReadSequence("CGTCCCCACAGACATCGGGT", ""), mappings);
//    for (auto &mapping: mappings) {
//        cout << mapping->to_string() << endl;
//    }

    logger->log_msg(LogLevel::INFO, "main", "Mapping reads...");
    read_mapping.do_mapping();

    auto reads = read_mapping.get_reads();
    long success = 0;
    for (auto &read: reads) {
        if (read->left->mapping != nullptr) {
            success++;
        }
        if (read->right->mapping != nullptr) {
            success++;
        }
    }
    logger->log_msg(LogLevel::INFO, "main", "Success: " + to_string(success)
                                            + "/" + to_string(reads.size() * 2) + " (" +
                                            to_string((double) success / (reads.size() * 2) * 100) + "%)");

    if (map_mode) {
        logger->log_msg(LogLevel::INFO, "main", "Outputting mappings...");
        read_mapping.output_mappings(output_file);
    } else {
        logger->log_msg(LogLevel::INFO, "main", "Calculating mutations...");
        read_mapping.calculate_mutations();

        logger->log_msg(LogLevel::INFO, "main", "Outputting mutations...");
        read_mapping.output_mutations(output_file);
    }
    logger->log_msg(LogLevel::INFO, "main", "Done.");

    return 0;
}
