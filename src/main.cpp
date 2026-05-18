#include "burnham/commands.hpp"

#include <iostream>

int main(int argc, char** argv) {
    return burnham::run_cli(argc, argv, std::cout, std::cerr);
}
