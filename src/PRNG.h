#pragma once
#include <random>
#include <chrono>

double urand(double lb, double up); // returns a uniformly distributed number between lb und up
void seed_nrand(unsigned int seed); // function to seed the XORSHIFT PRNG wit seed
void jump_nrand();                  // function to initialize a jump in the XORSHIFT PRNG
double nrand();                     // returns a normally distributed number


