#include <iostream>

#include "console.h"

void phonon_pumping::printFreqs(const list_of_frequencies& data) {
    for (std::size_t i = 0; i < data.size(); ++i) {
        for (const auto& [name, value] : data[i]) {
            std::cout << "  (" << name << ", " << value << " GHz)\n";
        }

        std::cout << '\n';
    }
}