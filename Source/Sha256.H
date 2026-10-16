#ifndef ABLSTEADY_SHA256_H_
#define ABLSTEADY_SHA256_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

namespace ablsteady {

//! SHA-256 (FIPS 180-4), for the checksums of reference archives. Host only.
class Sha256
{
public:
    Sha256 ();
    void update (const void* data, std::size_t n);
    //! The digest of everything given so far, as 64 lowercase hex
    //! characters. The object is unchanged and can take more data.
    [[nodiscard]] std::string hex () const;

private:
    std::string finish ();
    void block (const std::uint8_t* p);
    std::array<std::uint32_t, 8> m_h;
    std::array<std::uint8_t, 64> m_buf{};
    std::size_t m_buf_len{0};
    std::uint64_t m_total{0};
};

//! SHA-256 of a byte string, as hex.
std::string sha256_hex (const std::string& bytes);

//! SHA-256 of a file's contents, as hex. Returns an empty string if the file
//! cannot be opened or read.
std::string sha256_file (const std::string& path);

} // namespace ablsteady

#endif
