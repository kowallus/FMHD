#include "commons.hpp"

#include <iostream>

void
print_names(const std::vector<std::string> &names, const int len)
{
    *appout << "##Names ";
    for (int i = len - 1; i > -1;--i)
        *appout << names[i] << '\t';
    *appout << '\n';
}

std::ostream *appout = &std::cout;
FILE* appFile = stdout;
std::ostream *devout = &std::cout;

NullBuffer null_buffer;
std::ostream null_stream(&null_buffer);


// TIME

clock_t checkpoint;

clock_t clock_checkpoint() {
    checkpoint = clock();
    return checkpoint;
//    cout << "Clock reset!\n";
}

unsigned long long int clock_millis(clock_t checkpoint) {
    return (clock() - checkpoint) * (unsigned long long int) 1000 / CLOCKS_PER_SEC;
}

unsigned long long int clock_millis() {
    return clock_millis(checkpoint);
}


std::chrono::steady_clock::time_point chronocheckpoint;

std::chrono::steady_clock::time_point time_checkpoint() {
    chronocheckpoint = std::chrono::steady_clock::now();
    return chronocheckpoint;
}

unsigned long long int time_millis(std::chrono::steady_clock::time_point checkpoint) {
    std::chrono::nanoseconds time_span = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - checkpoint);
    return (double)time_span.count() / 1000000.0;
}

unsigned long long int time_millis() {
    return time_millis(chronocheckpoint);
}


struct LUT
{
    char complementsLut[256];
    char upperComplementsLut[256];

    char* cLutPtr = complementsLut - CHAR_MIN;
    char* uLutPtr = upperComplementsLut - CHAR_MIN;

    LUT() {
        for(int i = CHAR_MIN; i < CHAR_MAX; i++)
            uLutPtr[i] = i;
        uLutPtr['A'] = 'T'; uLutPtr['a'] = 'T';
        uLutPtr['C'] = 'G'; uLutPtr['c'] = 'G';
        uLutPtr['G'] = 'C'; uLutPtr['g'] = 'C';
        uLutPtr['T'] = 'A'; uLutPtr['t'] = 'A';
        uLutPtr['N'] = 'N'; uLutPtr['n'] = 'N';
        uLutPtr['U'] = 'A'; uLutPtr['u'] = 'A';
        uLutPtr['Y'] = 'R'; uLutPtr['y'] = 'R';
        uLutPtr['R'] = 'Y'; uLutPtr['r'] = 'Y';
        uLutPtr['K'] = 'M'; uLutPtr['k'] = 'M';
        uLutPtr['M'] = 'K'; uLutPtr['m'] = 'K';
        uLutPtr['B'] = 'V'; uLutPtr['b'] = 'V';
        uLutPtr['D'] = 'H'; uLutPtr['d'] = 'H';
        uLutPtr['H'] = 'D'; uLutPtr['h'] = 'D';
        uLutPtr['V'] = 'B'; uLutPtr['v'] = 'B';
        uLutPtr['W'] = 'S'; uLutPtr['w'] = 'S';
        uLutPtr['S'] = 'W'; uLutPtr['s'] = 'W';
        for(int i = CHAR_MIN; i < CHAR_MAX; i++)
            cLutPtr[i] = uLutPtr[i];
        cLutPtr['U'] = 'U';
        cLutPtr['u'] = 'u';
        cLutPtr['a'] = 't';
        cLutPtr['c'] = 'g';
        cLutPtr['g'] = 'c';
        cLutPtr['t'] = 'a';
        cLutPtr['n'] = 'n';
        cLutPtr['y'] = 'r';
        cLutPtr['r'] = 'y';
        cLutPtr['k'] = 'm';
        cLutPtr['m'] = 'k';
        cLutPtr['b'] = 'v';
        cLutPtr['d'] = 'h';
        cLutPtr['h'] = 'd';
        cLutPtr['v'] = 'b';
        cLutPtr['w'] = 's';
        cLutPtr['s'] = 'w';
    }

    void removeUpperWScomplements() {
        upperComplementsLut['W'] = 'W'; upperComplementsLut['w'] = 'w';
        upperComplementsLut['S'] = 'S'; upperComplementsLut['s'] = 's';
    }
} instance;

char* upperComplementsLUT = instance.uLutPtr;
char* complementsLUT = instance.cLutPtr;

char upperReverseComplement(char symbol) {
    return upperComplementsLUT[symbol];
}

char reverseComplement(char symbol) {
    return complementsLUT[symbol];
}

void upperReverseComplementInPlace(char* start, const std::size_t N) {
    char* left = start - 1;
    char* right = start + N;
    while (--right > ++left) {
        char tmp = upperComplementsLUT[*left];
        *left = upperComplementsLUT[*right];
        *right = tmp;
    }
    if (left == right)
        *left = upperComplementsLUT[*left];
}

void reverseComplementInPlace(char* start, const std::size_t N) {
    char* left = start - 1;
    char* right = start + N;
    while (--right > ++left) {
        char tmp = complementsLUT[*left];
        *left = complementsLUT[*right];
        *right = tmp;
    }
    if (left == right)
        *left = complementsLUT[*left];
}

void upperReverseComplement(const char* start, const std::size_t N, char* target) {
    char* right = target + N;
    while (right-- > target) {
        *right = upperComplementsLUT[*(start++)];
    }
}

void reverseComplement(const char* start, const std::size_t N, char* target) {
    char* right = target + N;
    while (right-- > target) {
        *right = complementsLUT[*(start++)];
    }
}

std::string upperReverseComplement(std::string kmer) {
    size_t kmer_length = kmer.size();
    std::string revcomp;
    revcomp.resize(kmer_length);
    size_t j = kmer_length;
    for(size_t i = 0; i < kmer_length; i++)
        revcomp[--j] = upperComplementsLUT[kmer[i]];
    return revcomp;
}

std::string reverseComplement(std::string kmer) {
    size_t kmer_length = kmer.size();
    std::string revcomp;
    revcomp.resize(kmer_length);
    size_t j = kmer_length;
    for(size_t i = 0; i < kmer_length; i++)
        revcomp[--j] = complementsLUT[kmer[i]];
    return revcomp;
}

void upperReverseComplementInPlace(std::string &kmer) {
    upperReverseComplementInPlace((char *) kmer.data(), kmer.length());
}

void reverseComplementInPlace(std::string &kmer) {
    reverseComplementInPlace((char *) kmer.data(), kmer.length());
}

void upperSequence(char* start, const std::size_t N) {
    char* left = start - 1;
    char* guard = start + N;
    while (++left < guard) {
        *left = toupper(*left);
    }
}
