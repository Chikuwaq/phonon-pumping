#include <iostream>

#include "console.h"

void phonon_pumping::printFreqs(const list_of_frequencies& data) {
    for (std::size_t i = 0; i < data.size(); ++i) {
        if (data[i].empty()) continue;
        
        for (const auto& [name, value] : data[i]) {
            std::cout << "  (" << name << ", " << value << " GHz)\n";
        }

        std::cout << '\n';
    }
}

void phonon_pumping::wait_for_enter() {
    std::cout << "Press Enter to continue..." << std::flush;
    std::string dummy;
    std::getline(std::cin, dummy);
}