#include "burnham/cigar.hpp"

#include <cctype>

namespace burnham {

Result<std::vector<CigarOp>> parse_cigar(const std::string& cigar) {
    std::vector<CigarOp> ops;
    if (cigar == "*") {
        return ops;
    }
    if (cigar.empty()) {
        return make_error("empty CIGAR");
    }

    std::int64_t length = 0;
    bool saw_digit = false;
    for (char ch : cigar) {
        if (std::isdigit(static_cast<unsigned char>(ch)) != 0) {
            saw_digit = true;
            length = (length * 10) + (ch - '0');
            continue;
        }
        const bool valid_op = ch == 'M' || ch == 'I' || ch == 'D' || ch == 'N' || ch == 'S' || ch == 'H' || ch == 'P' || ch == '=' || ch == 'X';
        if (!saw_digit || length <= 0 || !valid_op) {
            return make_error("invalid CIGAR: " + cigar);
        }
        ops.push_back(CigarOp{length, ch});
        length = 0;
        saw_digit = false;
    }
    if (saw_digit) {
        return make_error("trailing CIGAR length without op: " + cigar);
    }
    return ops;
}

CigarStats cigar_stats(const std::vector<CigarOp>& ops) {
    CigarStats stats;
    for (const auto& op : ops) {
        switch (op.op) {
        case 'M':
        case 'I':
        case 'S':
        case '=':
        case 'X':
            stats.query_length += op.length;
            break;
        default:
            break;
        }
        switch (op.op) {
        case 'M':
        case 'D':
        case 'N':
        case '=':
        case 'X':
            stats.reference_span += op.length;
            break;
        default:
            break;
        }
    }
    return stats;
}

} // namespace burnham
