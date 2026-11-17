#include <stdexcept>

#include "constants.h"

std::string phonon_pumping::get_mode_label(const enhance_mechanism mechanism) {
    std::string mode_label;
    if (mechanism == enhance_mechanism::EIGENMODE_TA) mode_label = "TA";
    else if (mechanism ==enhance_mechanism::EIGENMODE_LA) mode_label = "LA";
    else throw std::invalid_argument("Mode label cannot be inferred from 'enhance_mechanism' = " + std::to_string(mechanism));
    
    return mode_label;
}