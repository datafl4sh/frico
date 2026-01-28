/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2025, Matteo Cicuttin - IV3IWE
 * Politecnico di Torino
 * Dipartimento di Scienze Matematiche "G. L. Lagrange"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <expected>
#include <optional>
#include "utils.h"

namespace frico {

std::vector<std::string>
split(const std::string& str, char sep) {
    std::vector<std::string> tokens;

    size_t first = 0;
    size_t last = 0;

    while (first < str.length()) {

        while ( (first < str.length()) and str[first] == sep ) {
            first++;
        }

        last = first; 
        while ( (last < str.length()) and str[last] != sep ) {
            last++;
        }

        if (last > first) {
            tokens.push_back( str.substr(first, last - first) );
        }

        first = last+1;
    }

    return tokens;
}

std::expected<frequency_range, parse_error>
parse_frequency_range(const std::string& str)
{
    auto tokens = split(str, ':');
    if (tokens.size() != 3) {
        return std::unexpected(parse_error::invalid_input);
    }

    frequency_range ret;
    try {
        ret.start = std::stod(tokens[0]);
        ret.step = std::stod(tokens[1]);
        ret.end = std::stod(tokens[2]);
    }
    catch(...) {
        return std::unexpected(parse_error::invalid_input);
    }

    if ( (ret.start < 0.0) or (ret.step < 0.0) or (ret.end < 0.0) ) {
        return std::unexpected(parse_error::out_of_range);
    }

    if ( ret.start >= ret.end ) {
        return std::unexpected(parse_error::out_of_range);
    }

    return ret;
}

std::expected<std::vector<int>, parse_error>
parse_integer_list(const std::string& str)
{
    auto tokens = split(str, ',');

    std::vector<int> ret;
    ret.reserve( tokens.size() );

    for (auto& tok : tokens) {
        try {
            ret.push_back( std::stoi(tok) );
        }
        catch (...) {
            return std::unexpected(parse_error::invalid_input);
        }
    }

    return ret;
}

std::optional<frico::frequency_range>
parse_frequency_parameters(const char *arg_frequency, const char *arg_range_expr)
{
    if ( (arg_frequency == nullptr) and (arg_range_expr == nullptr) ) {
        std::cerr << "Simulation frequency not specified (-f or -R)\n";
        return {};
    }

    if ( (arg_frequency != nullptr) and (arg_range_expr != nullptr) ) {
        std::cerr << "Flags -f and -R are mutually exclusive\n";
        return {};
    }

    if ( arg_frequency ) {
        try {
            double frequency = std::stod(arg_frequency);
            return frico::frequency_range {
                .start = frequency,
                .step = frequency,
                .end = frequency
            };
        }
        catch (...) {
            std::cerr << "Error parsing the argument of -f\n";
            return {};
        }
    }

    if ( arg_range_expr ) {
        auto exp_range = frico::parse_frequency_range(arg_range_expr);
        if ( exp_range.has_value() ) {
            return *exp_range;
        }
        
        switch ( exp_range.error() ) {
            case frico::parse_error::invalid_input:
                std::cerr << "Malformed range expression for -R\n";
                return {};
            case frico::parse_error::out_of_range:
                std::cerr << "Invalid range for -R\n";
                return {};
            default:
                return {};
        }
    }

    return {};
}

} // namespace frico