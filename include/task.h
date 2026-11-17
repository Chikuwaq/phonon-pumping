/*
 * SPDX-FileCopyrightText: 2025 Takuma Sato
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include "constants.h"

namespace phonon_pumping {

    namespace task {
        enum class TaskID {
            PLOT_FMR_SPECTRUM,
            FIT_TA_SOUND_VELOCITIES,
            FIT_LA_SOUND_VELOCITIES
        };

        std::vector<TaskID> load(const std::filesystem::path& filepath);
        std::vector<double> extract_freqs(const std::filesystem::path& filepath, const enhance_mechanism mode);
    }

}