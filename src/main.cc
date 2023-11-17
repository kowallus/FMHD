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
    const int max_reading_threads = 4;
    int num_of_threads;
#pragma omp parallel
#pragma omp single
    num_of_threads = omp_get_num_threads();
    app.add_option("-c,--cmd", cmd, "mode in which program operates")->required();
    app.add_option("-l,--inputlist", files_list, "filenames list for comparison");
    app.add_option("-k,--kmerlen", kmerlen, "k-mer lenght. Default: 21");
    app.add_option("-o,--outfile", out_file, "output filename with sketches. Default: sketch.bin");
    app.add_option("-t,--nthreads", num_of_threads, "number of threads used for comparison.");
    app.add_flag("-E,--edgelist", edgelist, "return output as a edge list with fields [seq1, seq2, dist], default: triangle matrix");
    CLI11_PARSE(app, argc, argv);

    omp_set_dynamic(false);
    std::chrono::steady_clock::time_point start_t = std::chrono::steady_clock::now();

//    appout = &null_stream;

    if (cmd == "sketch") {
        std::string line;
        std::ifstream input_file(files_list);
        std::vector<std::string> files_names;
        std::vector<std::string> genomes;

        while (std::getline(input_file, line)) {
            files_names.push_back(line);
        }
        const size_t len_sketches = files_names.size();
        genomes.resize(len_sketches);
        int reading_threads = num_of_threads > max_reading_threads ? max_reading_threads : num_of_threads;
        omp_set_num_threads(reading_threads);
        *devout << "number of threads: " << reading_threads << std::endl;
        #pragma omp parallel for
        for (int i = 0; i < len_sketches; ++i) {
            genomes[i] = read_fna(files_names[i]);
        }
        *devout << "reading in "
                << time_millis(start_t) << " msec" << std::endl;
        start_t = std::chrono::steady_clock::now();

        omp_set_num_threads(num_of_threads);
        *devout << "number of threads: " << num_of_threads << std::endl;
        std::vector<std::array<uint64_t, M>> sketches(len_sketches);

        #pragma omp parallel for
        for (int i = 0; i < len_sketches; ++i) {
            sketches[i] = get_sketch_hash_template<MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>(genomes[i], kmerlen);
//            sketches[i] = get_sketch(genomes[i], (uint8_t) kmerlen, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID); // slower
//            sketches[i] = get_sketch_full_template<8, NONE_HASH_ID>(genomes[i]); // only for testing, slower than without templates
        }

        save_sketches(sketches, files_names, out_file);
    } else {
        omp_set_num_threads(num_of_threads);
        *devout << "number of threads: " << num_of_threads << std::endl;

        std::vector<std::array<uint64_t, M>> sketches(0);
        std::vector<std::string> files_names(0);
        if (read_sketches(sketches, files_names, files_list))
            functionMap.at(cmd)(sketches, files_names, num_of_threads, edgelist);
    }

    *devout << "finished in "
            << time_millis(start_t) << " msec" << std::endl;

    return EXIT_SUCCESS;
}
