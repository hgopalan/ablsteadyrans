// Beyond the plan: SHA-256 against the FIPS 180-4 / NIST example vectors.
// The reference archives (UT-1.1.1) rest on it.
#include "Sha256.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

namespace {

TEST (Sha256, NistVectors)
{
    EXPECT_EQ(ablsteady::sha256_hex(""),
              "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    EXPECT_EQ(ablsteady::sha256_hex("abc"),
              "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    EXPECT_EQ(ablsteady::sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
              "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    EXPECT_EQ(ablsteady::sha256_hex(std::string(1000000, 'a')),
              "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

// Padding boundaries: messages of 55, 56, 63, 64 and 65 bytes, fed whole and
// one byte at a time, give the same digest.
TEST (Sha256, ChunkingDoesNotMatter)
{
    for (const std::size_t n : {55u, 56u, 63u, 64u, 65u, 1000u}) {
        std::string msg(n, '\0');
        for (std::size_t i = 0; i < n; ++i) { msg[i] = static_cast<char>((i * 31 + 7) & 0xff); }
        ablsteady::Sha256 s;
        for (const char c : msg) { s.update(&c, 1); }
        EXPECT_EQ(s.hex(), ablsteady::sha256_hex(msg)) << n;
    }
}

// hex() leaves the object usable: asking twice gives the same digest, and
// more data afterwards gives the digest of the longer message.
TEST (Sha256, HexIsRepeatable)
{
    ablsteady::Sha256 s;
    s.update("ab", 2);
    EXPECT_EQ(s.hex(), s.hex());
    s.update("c", 1);
    EXPECT_EQ(s.hex(), ablsteady::sha256_hex("abc"));
}

TEST (Sha256, FileDigestMatchesBytes)
{
    const std::string name = "ut_sha256_" + std::to_string(amrex::ParallelDescriptor::MyProc()) + ".bin";
    std::string msg(200000, '\0');
    for (std::size_t i = 0; i < msg.size(); ++i) { msg[i] = static_cast<char>(i * 2654435761u >> 24); }
    {
        std::ofstream os(name, std::ios::binary);
        os.write(msg.data(), static_cast<std::streamsize>(msg.size()));
    }
    EXPECT_EQ(ablsteady::sha256_file(name), ablsteady::sha256_hex(msg));
    std::remove(name.c_str());
    EXPECT_EQ(ablsteady::sha256_file(name), "");   // missing file
}

} // namespace
