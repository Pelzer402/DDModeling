#include "PRNG.h"

uint64_t s[2];

uint64_t seed_seeder;
std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());

double urand(double lb, double up)
{
  std::uniform_real_distribution<double> distribution(lb,up);
	return distribution(generator);
}
uint64_t Seeder_xor() {
	uint64_t z = (seed_seeder += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

void seed_nrand(unsigned int seed)
{
	seed_seeder = seed;
	s[0] = Seeder_xor();
	s[1] = Seeder_xor();
}
static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}
double next(void) {
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;
	s1 ^= s0;
	s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
	s[1] = rotl(s1, 36); // c
	return (result >> 11) * (1. / (UINT64_C(1) << 53));
}
void jump_nrand()
{
	static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for (unsigned long long int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for (int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			next();
		}
	s[0] = s0;
	s[1] = s1;
}

double nrand()
{ // Generating unit normal distributed random variables by the polar method
	double u1, u2, v1, v2, s, z;
	static double gset;
	static int iset = 0;

	if (iset == 0) {
		do
		{
			u1 = next();
			u2 = next();
			v1 = 2.0*u1 - 1;
			v2 = 2.0*u2 - 1;
			s = v1*v1 + v2*v2;
		} while (s >= 1 || s == 0);
		z = sqrt(-2 * log(s) / s);
		gset = z*v1;
		iset = 1;
		return z*v2;
	}
	else {
		iset = 0;
		return gset;
	}
}

