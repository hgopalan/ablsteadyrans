#ifndef GUST_FILE_MAGIC_H_
#define GUST_FILE_MAGIC_H_

#include <string>

namespace gust {

//! True if `word`, the first word of a state file, is `magic` ("gust-<kind>")
//! or the "ablsteady-<kind>" written before the solver was renamed GUST
//! (October 2026), so older checkpoints and histories still restart.
inline bool is_file_magic (const std::string& word, const std::string& magic)
{
    const std::string prefix = "gust-";
    return word == magic
        || (magic.compare(0, prefix.size(), prefix) == 0
            && word == "ablsteady-" + magic.substr(prefix.size()));
}

} // namespace gust

#endif
