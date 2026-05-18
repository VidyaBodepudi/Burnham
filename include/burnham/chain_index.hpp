#pragma once

#include "burnham/result.hpp"

#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace burnham {

enum class SourceSide { Query, Target };

struct ReferenceDictionary {
    std::map<std::string, std::int64_t> lengths;
};

struct ChainBlock {
    std::string source_contig;
    std::int64_t source_start = 0;
    std::int64_t source_end = 0;
    std::string dest_contig;
    std::int64_t dest_start = 0;
    std::int64_t dest_end = 0;
    char strand = '+';
    std::string chain_id;
};

struct ChainIndex {
    SourceSide source_side = SourceSide::Query;
    std::vector<ChainBlock> blocks;
};

struct MappedPosition {
    std::string contig;
    std::int64_t position = 0;
    char strand = '+';
    std::string chain_id;
};

struct MappedInterval {
    std::string contig;
    std::int64_t start = 0;
    std::int64_t end = 0;
    char strand = '+';
    std::string chain_id;
};

Result<ReferenceDictionary> read_fai(const std::filesystem::path& path);
Result<ChainIndex> build_chain_index(const std::filesystem::path& chain_path,
                                     SourceSide source_side,
                                     const ReferenceDictionary* source_dict,
                                     const ReferenceDictionary* dest_dict);
Result<void> write_chain_index(const ChainIndex& index, const std::filesystem::path& path, bool force);
Result<ChainIndex> read_chain_index(const std::filesystem::path& path);
std::optional<MappedPosition> map_position(const ChainIndex& index, const std::string& contig, std::int64_t position);
std::optional<MappedInterval> map_interval(const ChainIndex& index, const std::string& contig, std::int64_t start, std::int64_t end);
std::string inspect_chain_index_text(const ChainIndex& index, bool as_json);
std::string source_side_name(SourceSide side);
Result<SourceSide> parse_source_side(const std::string& value);

} // namespace burnham
