#include <iostream> /* display output in console */
#include <array> /* array type used in sketches, mixers, etc */
#include <omp.h>

#include "commons.hpp"
#include "constants.hpp"
#include "sketch.hpp"
#include "CLI11.hpp" /* parsing command-line arguments */
#include "sketch_file_handler.hpp"
#include "distances/dist_map.hpp"

std::string read_fna(const std::string& filename)
{
    std::ifstream ifs (filename, std::ifstream::in);
    std::string line="", output = "";
    while (std::getline(ifs, line)) {
        if (line.empty() || line.at(0) == '>')
            continue;
        output.append(line);
    }
    ifs.close();
    return output;
}

int main (int argc, char **argv) {
    CLI::App app{"library for comparing genomes"};
    std::string files_list;
    std::string out_file = "sketch.bin";
    std::string cmd = "dist";
    size_t kmerlen = 21;
    bool edgelist = false;
    int num_of_threads = 1;
    app.add_option("-c,--cmd", cmd, "mode in which program operates")->required();
    app.add_option("-l,--inputlist", files_list, "filenames list for comparison");
    app.add_option("-k,--kmerlen", kmerlen, "k-mer lenght. Default: 21");
    app.add_option("-o,--outfile", out_file, "output filename with sketches. Default: sketch.bin");
    app.add_option("-t,--nthreads", num_of_threads, "number of threads used for comparison. Default: 1");
    app.add_flag("-E,--edgelist", edgelist, "return output as a edge list with fields [seq1, seq2, dist], default: triangle matrix");
    CLI11_PARSE(app, argc, argv);

    omp_set_dynamic(false);
    omp_set_num_threads(num_of_threads);

    std::chrono::steady_clock::time_point start_t = std::chrono::steady_clock::now();

//    appout = &null_stream;

    if (cmd == "sketch") {
        std::string line;
        std::ifstream input_file(files_list);
        std::vector<std::string> files_names;
        std::vector<std::string> genomes;

        while (std::getline(input_file, line)) {
            files_names.push_back(line);
            genomes.push_back(read_fna(line));
        }
        const size_t len_sketches = files_names.size();

        *devout << "reading in "
                << time_millis(start_t) << " msec" << std::endl;
        start_t = std::chrono::steady_clock::now();

        std::vector<std::array<uint64_t, M>> sketches(len_sketches);

        #pragma omp parallel for
        for (int i = 0; i < len_sketches; ++i) {
            sketches[i] = get_sketch_hash_template<MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>(genomes[i], kmerlen);
//            sketches[i] = get_sketch(genomes[i], (uint8_t) kmerlen, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID); // slower
//            sketches[i] = get_sketch_full_template<4, NONE_HASH_ID>(genomes[i]); // only for testing, slower than without templates
        }

        save_sketches(sketches, files_names, out_file);
    } else {
        std::vector<std::array<uint64_t, M>> sketches(0);
        std::vector<std::string> files_names(0);
        if (read_sketches(sketches, files_names, files_list))
            functionMap.at(cmd)(sketches, files_names, num_of_threads, edgelist);
    }

    *devout << "finished in "
            << time_millis(start_t) << " msec" << std::endl;

    return EXIT_SUCCESS;
}
