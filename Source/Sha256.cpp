#include "Sha256.H"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <vector>

namespace ablsteady {

namespace {

constexpr std::array<std::uint32_t, 64> kK{
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};

constexpr std::uint32_t rotr (std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

} // namespace

Sha256::Sha256 ()
    : m_h{0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
          0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u}
{}

void Sha256::block (const std::uint8_t* p)
{
    std::array<std::uint32_t, 64> w{};
    for (int i = 0; i < 16; ++i) {
        w[i] = (std::uint32_t(p[4 * i]) << 24) | (std::uint32_t(p[4 * i + 1]) << 16)
             | (std::uint32_t(p[4 * i + 2]) << 8) | std::uint32_t(p[4 * i + 3]);
    }
    for (int i = 16; i < 64; ++i) {
        const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    std::uint32_t a = m_h[0], b = m_h[1], c = m_h[2], d = m_h[3];
    std::uint32_t e = m_h[4], f = m_h[5], g = m_h[6], h = m_h[7];
    for (int i = 0; i < 64; ++i) {
        const std::uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        const std::uint32_t ch = (e & f) ^ (~e & g);
        const std::uint32_t t1 = h + S1 + ch + kK[i] + w[i];
        const std::uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        const std::uint32_t t2 = S0 + maj;
        h = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    m_h[0] += a; m_h[1] += b; m_h[2] += c; m_h[3] += d;
    m_h[4] += e; m_h[5] += f; m_h[6] += g; m_h[7] += h;
}

void Sha256::update (const void* data, std::size_t n)
{
    const auto* p = static_cast<const std::uint8_t*>(data);
    m_total += n;
    while (n > 0) {
        const std::size_t take = std::min<std::size_t>(64 - m_buf_len, n);
        for (std::size_t i = 0; i < take; ++i) { m_buf[m_buf_len + i] = p[i]; }
        m_buf_len += take;
        p += take;
        n -= take;
        if (m_buf_len == 64) {
            block(m_buf.data());
            m_buf_len = 0;
        }
    }
}

std::string Sha256::hex () const
{
    Sha256 copy(*this);
    return copy.finish();
}

std::string Sha256::finish ()
{
    const std::uint64_t bits = m_total * 8;
    const std::uint8_t one = 0x80;
    update(&one, 1);
    const std::uint8_t zero = 0;
    while (m_buf_len != 56) { update(&zero, 1); }
    std::uint8_t len[8];
    for (int i = 0; i < 8; ++i) { len[i] = static_cast<std::uint8_t>(bits >> (56 - 8 * i)); }
    update(len, 8);

    std::string out;
    char buf[9];
    for (const auto v : m_h) {
        std::snprintf(buf, sizeof(buf), "%08x", v);
        out += buf;
    }
    return out;
}

std::string sha256_hex (const std::string& bytes)
{
    Sha256 s;
    s.update(bytes.data(), bytes.size());
    return s.hex();
}

std::string sha256_file (const std::string& path)
{
    std::ifstream is(path, std::ios::binary);
    if (!is) { return {}; }
    Sha256 s;
    std::vector<char> buf(1 << 16);
    while (is) {
        is.read(buf.data(), static_cast<std::streamsize>(buf.size()));
        const auto got = is.gcount();
        if (got > 0) { s.update(buf.data(), static_cast<std::size_t>(got)); }
    }
    if (!is.eof()) { return {}; }
    return s.hex();
}

} // namespace ablsteady
