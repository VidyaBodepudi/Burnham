#pragma once

#include "burnham/result.hpp"

#include <iosfwd>

namespace burnham {

int run_cli(int argc, char** argv, std::ostream& out, std::ostream& err);

} // namespace burnham
