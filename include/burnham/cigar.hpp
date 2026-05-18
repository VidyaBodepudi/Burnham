#pragma once

#include "burnham/result.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace burnham {

struct CigarOp {
    std::int64_t length = 0;
    char op = 'M';
};

struct CigarStats {
    std::int64_t query_length = 0;
    std::int64_t reference_span = 0;
};

Result<std::vector<CigarOp>> parse_cigar(const std::string& cigar);
CigarStats cigar_stats(const std::vector<CigarOp>& ops);

} // namespace burnham
