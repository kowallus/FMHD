#pragma once

#ifndef _COMMONS_H_
#define _COMMONS_H_

#include <string>
#include <vector>
#include <streambuf>
#include <chrono>
#include <climits>

void print_names(const std::vector<std::string> &names, const int len);

class NullBuffer : public std::streambuf
{
public:
    int overflow(int c) { return c; }
};

extern std::ostream null_stream;

// time routines

clock_t clock_checkpoint();
unsigned long long int clock_millis();
unsigned long long int clock_millis(clock_t checkpoint);

std::chrono::steady_clock::time_point time_checkpoint();
unsigned long long int time_millis();
unsigned long long int time_millis(std::chrono::steady_clock::time_point checkpoint);

// input output routines

extern std::ostream *appout;
extern std::ostream *devout;

// bioinformatical routines

char upperReverseComplement(char symbol);
char reverseComplement(char symbol);
std::string upperReverseComplement(std::string kmer);
std::string reverseComplement(std::string kmer);
void upperReverseComplement(const char* start, const std::size_t N, char* target);
void reverseComplement(const char* start, const std::size_t N, char* target);
void upperReverseComplementInPlace(char* start, const std::size_t N);
void reverseComplementInPlace(char* start, const std::size_t N);
void upperReverseComplementInPlace(std::string &kmer);
void reverseComplementInPlace(std::string &kmer);
void upperSequence(char* start, const std::size_t N);

#endif