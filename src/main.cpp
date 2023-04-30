#include <iostream>

#include "ReadMapping.h"
#include "logging.h"

using namespace std;

int main(int argc, char *argv[]) {
    // Parse command line arguments
    string reference_file, input_file;
    if (argc == 3) {
        reference_file = argv[1];
        input_file = argv[2];
    } else {
        cout << "Usage: ./release reference_file input_file" << endl;
        return 1;
    }

    auto logger = new Logger(LogLevel::INFO_ONLY);

    logger->log_msg(LogLevel::INFO, "main", "Reading input files...");
    ReadMapping read_mapping(reference_file, input_file, logger);

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
                                            + "/" + to_string(reads.size()) + " (" +
                                            to_string((double) success / (reads.size() * 2) * 100) + "%)");

    logger->log_msg(LogLevel::INFO, "main", "Calculating mutations...");
    read_mapping.calculate_mutations();

    logger->log_msg(LogLevel::INFO, "main", "Outputting mutations...");
    read_mapping.output_mutations("output/output.txt");

    logger->log_msg(LogLevel::INFO, "main", "Done.");

    return 0;
}
