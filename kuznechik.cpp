#pragma GCC optimize ("unroll-loops")
#include <iostream>
#include <iomanip>
#include <vector>
#include <bitset>
#include <stdint.h>
#include <charconv>
#include <string>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <chrono>
#include <random>
#include <array>


namespace NField {
    const uint16_t MOD = 0b111000011;
    std::vector<uint8_t> log(256);
    std::vector<uint8_t> exp(256);
    uint16_t SlowMultiply(uint16_t first, uint16_t second) {
        uint16_t result = 0;
        for (int i = 0; i < 8; ++i) {
            if (first & (1 << i)) {
                result ^= second * (1 << i);
            }
        }
        for (int i = 14; i >= 8; --i) {
            if (result & (1 << i)) {
                result ^= MOD * (1 << (i - 8));
            }
        }
        return result;
    }
    void Init() {
        log.resize(256);
        exp.resize(256);
        uint16_t root = 2;
        uint16_t current = 1;
        log[current] = 0;
        exp[0] = current;
        for (int i = 1; i < 255; ++i) {
            current = SlowMultiply(current, root);
            log[current] = i;
            exp[i] = current;
        }
    }
    class TElement {
        uint8_t body;
    public:
        TElement() : body(0) {}
        explicit TElement(uint8_t value) : body(value) {}
        explicit TElement(int value) : body(value) {}
        uint8_t GetBody() const {
            return this->body;
        }
        void SetBody(uint8_t value) {
            this->body = value;
        }
        TElement& operator*=(const TElement& other) {
            if (this->body == 0 || other.body == 0) {
                this->body = 0;
                return *this;
            }
            this->body = exp[(log[this->body] + log[other.body]) % 255];
            return *this;
        }
        TElement& operator+=(const TElement& other) {
            this->body ^= other.body;
            return *this;
        }
        TElement& operator^=(const uint8_t& value) {
            this->body ^= value;
            return *this;
        }
    };

    TElement operator +(TElement left, const TElement& right) {
        left += right;
        return left;
    }

    TElement operator *(TElement left, const TElement& right) {
        left *= right;
        return left;
    }
}


namespace NKuznechik {
    typedef NField::TElement TBlock[16];
    typedef uint8_t TBytes[16];
    // TBlock Bytes2Block(const TBytes& bytes) {
    //     if (bytes.size() != 16) {
    //         throw std::invalid_argument("Bad size");
    //     }
    //     TBlock block;
    //     for (int i = 0; i < 16; ++i) {   
    //         block[i] = NField::TElement(bytes[i]);
    //     }
    //     return block;
    // }
    const NField::TElement GAMMA_SEQUENCE[16] = {
        NField::TElement(148), NField::TElement(32), NField::TElement(133), NField::TElement(16),
        NField::TElement(194), NField::TElement(192), NField::TElement(1), NField::TElement(251),
        NField::TElement(1), NField::TElement(192), NField::TElement(194), NField::TElement(16),
        NField::TElement(133), NField::TElement(32), NField::TElement(148), NField::TElement(1)
    };
    const uint8_t PERMUTATION[256] = {
        252,238,221,17,207,110,49,22,251,196,250,218,35,197,4,77,
        233,119,240,219,147,46,153,186,23,54,241,187,20,205,95,193,
        249,24,101, 90,226,92,239,33,129,28,60,66,139,1,142,79,
        5,132,2,174,227,106,143, 160,6,11,237,152,127,212,211,31,
        235,52,44,81,234,200,72,171,242,42,104,162,253,58,206,204,
        181,112,14,86,8,12,118,18,191,114,19,71,156,183,93,135,
        21,161,150,41,16,123,154,199,243,145,120,111,157,158,178,177,
        50,117,25,61,255,53,138,126,109,84,198,128,195,189,13,87,
        223,245,36,169,62,168,67,201,215,121,214,246,124,34,185,3,
        224,15,236,222,122,148,176,188,220,232,40,80,78,51,10,74,
        167,151,96,115,30,0,98,68,26,184,56,130,100,159,38,65,
        173,69,70,146,39,94,85,47,140,163,165,125,105,213,149,59,
        7,88,179,64,134,172,29,247,48,55,107,228,136,217,231,137,
        225,27,131,73,76,63,248,254,141,83,170,144,202,216,133,97,
        32,113,103,164,45,43,9,91,203,155,37,208,190,229,108,82,
        89,166,116,210,230,244,180,192,209,102,175,194,57,75,99,182
    };
    uint16_t S[65536];
    uint16_t REVERSED_S[65536];
    uint8_t ITERATION_CONSTANTS[32 * 16];
    NField::TElement LINEAR_PRECALC[16 * 256 * 16];
    NField::TElement REVERSED_LINEAR_PRECALC[16 * 256 * 16];
    NField::TElement Gamma(NField::TElement* block, int block_start) {
        NField::TElement result;
        for (int i = 0; i < 16; ++i) {
            result += block[(block_start + i) % 16] * GAMMA_SEQUENCE[15 - i];
        }
        return result;
    }
    void SlowH(NField::TElement* block, bool inversed = false) {
        int block_start = 0;
        for (int i = 0; i < 16; ++i) {
            NField::TElement extra;
            if (inversed) {
                block_start = (block_start - 1 + 16) % 16;
                block[block_start] = Gamma(block, block_start);
            }
            else {
                block[block_start] = Gamma(block, block_start);
                block_start = (block_start + 1) % 16;
            }
        }
        NField::TElement tmp[16];
        for (int i = 0; i < 16; ++i) {
            tmp[i] = block[(i + block_start) % 16];
        }
        std::memcpy(block, tmp, 16);
    }
    void Add(NField::TElement* block, uint8_t* key) {
        uint64_t* roflan1 = (uint64_t*)block;
        uint64_t* roflan2 = (uint64_t*)block + 1;

        uint64_t* other1 = (uint64_t*)key;
        uint64_t* other2 = (uint64_t*)key + 1;

        *roflan1 ^= *other1;
        *roflan2 ^= *other2;
    }
    void Add16(NField::TElement* block, NField::TElement* other) {

        uint64_t* roflan1 = (uint64_t*)block;
        uint64_t* roflan2 = (uint64_t*)block + 1;

        uint64_t* other1 = (uint64_t*)other;
        uint64_t* other2 = (uint64_t*)other + 1;

        *roflan1 ^= *other1;
        *roflan2 ^= *other2;
    }
    void H(NField::TElement* block, bool inversed = false) {
        NField::TElement result[16];
        if (!inversed) {
            for (int i = 0; i < 16; ++i) {
                Add16(result, LINEAR_PRECALC + (256 * i + block[i].GetBody()) * 16);
            }
        }
        else {
            for (int i = 0; i < 16; ++i) {
                Add16(result, REVERSED_LINEAR_PRECALC + (256 * i + block[i].GetBody()) * 16);
            }
        }
        std::memcpy(block, result, 16);
    }
    void N(NField::TElement* block, bool inversed = false) {
        uint16_t* ptr = (uint16_t*)block;
        if (inversed) {
            (*ptr) = REVERSED_S[(*ptr)];
            (*(ptr + 1)) = REVERSED_S[(*(ptr + 1))];
            (*(ptr + 2)) = REVERSED_S[(*(ptr + 2))];
            (*(ptr + 3)) = REVERSED_S[(*(ptr + 3))];
            (*(ptr + 4)) = REVERSED_S[(*(ptr + 4))];
            (*(ptr + 5)) = REVERSED_S[(*(ptr + 5))];
            (*(ptr + 6)) = REVERSED_S[(*(ptr + 6))];
            (*(ptr + 7)) = REVERSED_S[(*(ptr + 7))];
        }
        else {
            (*ptr) = S[(*ptr)];
            (*(ptr + 1)) = S[(*(ptr + 1))];
            (*(ptr + 2)) = S[(*(ptr + 2))];
            (*(ptr + 3)) = S[(*(ptr + 3))];
            (*(ptr + 4)) = S[(*(ptr + 4))];
            (*(ptr + 5)) = S[(*(ptr + 5))];
            (*(ptr + 6)) = S[(*(ptr + 6))];
            (*(ptr + 7)) = S[(*(ptr + 7))];
        }
    }
    void KeyGeneration(TBytes* storage, uint8_t* key) {
        NField::TElement k1[16], k2[16];
        std::memcpy(k2, key, 16);
        std::memcpy(k1, key + 16, 16);
        std::memcpy(storage, k1, 16);
        std::memcpy(storage + 1, k2, 16);
        NField::TElement tmp[16];
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 8; ++i) {
                std::memcpy(tmp, k1, 16);
                Add(k1, ITERATION_CONSTANTS + (8 * j + i) * 16);
                N(k1);
                H(k1);
                Add16(k1, k2);
                std::memcpy(k2, tmp, 16);
            }
            std::memcpy(storage + 2 * (j + 1), k1, 16);
            std::memcpy(storage + 2 * (j + 1) + 1, k2, 16);
        }
    }
    TBytes kkeys[10];
    bool are = false;
    void Encode(NField::TElement* block, uint8_t* key) {
        if (!are) {
            KeyGeneration(kkeys, key);
            are = true;
        }
        for (int i = 0; i < 9; ++i) {
            Add(block, kkeys[i]);
            N(block);
            H(block);
        }
        Add(block, kkeys[9]);
    }
    void Decode(NField::TElement* block, uint8_t* key) {
        if (!are) {
            KeyGeneration(kkeys, key);
            are = true;
        }
        for (int i = 9; i > 0; --i) {
            Add(block, kkeys[i]);
            H(block, true);
            N(block, true);
        }
        Add(block, kkeys[0]);
    }
    void PrintBlock(NField::TElement* block) {
        for (int i = 0; i < 16; ++i) {
            std::cout << std::hex << std::setw(2) << std::setfill('0') << std::right << int(block[15 - i].GetBody());
        }
        std::cout << "\n";
    }
    void KeyString2Bytes(uint8_t* storage, const std::string& s) {
        if (s.length() != 64) {
            throw std::invalid_argument("Bad size");
        }
        int n = 32;
        for (int i = 0; i + 1 < 64; i += 2) {
            storage[n - 1 - i / 2] = std::stoi(s.substr(i, 2), nullptr, 16);
        }
    }
    void String2Block(NField::TElement* storage, const std::string& s) {
        if (s.length() != 32) {
            throw std::invalid_argument("Bad size");
        }
        for (int i = 0; i < 16; ++i) {
            storage[15 - i].SetBody(std::stoi(s.substr(2 * i, 2), nullptr, 16));
        }
    }
    void Init() {
        for (int i = 0; i < 65536; ++i) {
            uint16_t from = i, to = 0;
            uint8_t* ptr_from = (uint8_t*)(&from);
            uint8_t* ptr_to = (uint8_t*)(&to);
            (*ptr_to) = PERMUTATION[*ptr_from];
            (*(ptr_to + 1)) = PERMUTATION[*(ptr_from + 1)];
            S[from] = to;
            REVERSED_S[to] = from;
        }
        for (int i = 0; i < 32; ++i) {
            NField::TElement a[16];
            a[0].SetBody(i + 1);
            SlowH(a);
            for (int j = 0; j < 16; ++j) {
                ITERATION_CONSTANTS[i * 16 + j] = a[j].GetBody();
            }
        }
        for (int i = 0; i < 16; ++i) {
            for (int j = 0; j < 256; ++j) {
                {
                    NField::TElement block[16];
                    block[i].SetBody(j);
                    SlowH(block);
                    std::memcpy(LINEAR_PRECALC + (256 * i + j) * 16, block, 16);
                }
                {
                    NField::TElement block[16];
                    block[i].SetBody(j);
                    SlowH(block, true);
                    std::memcpy(REVERSED_LINEAR_PRECALC + (256 * i + j) * 16, block, 16);
                }
            }
        }
    }
}


int main() {
    std::ios::sync_with_stdio(false); std::cin.tie(0);
    freopen("in.txt", "r", stdin);
    /*for (int i = 0; i < 1000000; ++i) {
        for (int j = 0; j < 16; ++j) {
            unsigned int x = rand() % 256;
            std::cout << std::hex << std::setw(2) << std::setfill('0') << std::right << x;
        }
        std::cout << "\n";
    }*/
    NField::Init();
    NKuznechik::Init();

    uint8_t key[32];
    NKuznechik::KeyString2Bytes(key, "8899aabbccddeeff0011223344556677fedcba98765432100123456789abcdef");

    // Encoding speed test
    {
        int blocks_to_read = 1000000;
        std::vector<NField::TElement*> s(blocks_to_read);
        int cnt = 0;
        for (int i = 0; i < blocks_to_read; ++i) {
            std::string ss;
            std::cin >> ss;
            s[i] = new NField::TElement[16];
            NKuznechik::String2Block(s[i], ss);
        }
        auto begin = std::chrono::high_resolution_clock::now();
        NField::TElement block[16];
        for (int i = 0; i < blocks_to_read; ++i) {
            ++cnt;
            NKuznechik::Encode(s[i], key);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = (long double)std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000;
        std::cout << "Time in seconds: " << elapsed << std::endl;
        std::cout << "Speed in MB/s: " << (long double)cnt * 16 / 1024 / 1024 / elapsed << std::endl;
        for (int i = 0; i < blocks_to_read; ++i) {
            delete[] s[i];
        }
    }

    // To test correctness
    {
        NField::TElement block[16];
        NKuznechik::String2Block(block, "1122334455667700ffeeddccbbaa9988");
        NKuznechik::Encode(block, key);
        std::ostringstream sout;
        for (int i = 0; i < 16; ++i) {
            sout << std::hex << std::setw(2) << std::setfill('0') << std::right << int(block[15 - i].GetBody());
        }
        if (sout.str() != "7f679d90bebc24305a468d42b9d4edcd") {
            std::cout << "Bad encoding" << std::endl;
        }
        NKuznechik::String2Block(block, "7f679d90bebc24305a468d42b9d4edcd");
        NKuznechik::Decode(block, key);
        sout.str("");
        sout.clear();
        for (int i = 0; i < 16; ++i) {
            sout << std::hex << std::setw(2) << std::setfill('0') << std::right << int(block[15 - i].GetBody());
        }
        if (sout.str() != "1122334455667700ffeeddccbbaa9988") {
            std::cout << "Bad decoding" << std::endl;
        }
    }
    return 0;
}
