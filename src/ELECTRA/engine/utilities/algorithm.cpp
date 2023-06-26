/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */


#include "ELECTRA/engine/utilities/algorithm.hpp"


namespace ELECTRA {

namespace ALGORITHM {


bool StringCompCaseInsensitive(const std::string &str1, const std::string &str2)
{
    if (str1.size() != str2.size())
        return false;

    return std::equal(std::begin(str1), std::end(str1), std::begin(str1), std::end(str1),
                      [](const char &a, const char &b) { return std::tolower(a) == std::tolower(b); });
}


bool ExistExactWord(const std::string &phrase, const std::string &word)
{
    // Guess it is found.
    bool found = true;

    // First character position.
    std::size_t cp_pos = phrase.find(word[0]);
    if(cp_pos == std::string::npos) { return false; }

    // Next characters position.
    std::size_t cn_pos = std::string::npos;
    for (std::size_t i = 1; i != word.size(); ++i) {

        // Search for the charecter starting from the next char than the previously found.
        cn_pos = phrase.find(word[i], cp_pos+1);

        // If not found or not next to the previous one, return false.
        if(cn_pos == std::string::npos || cn_pos != (cp_pos+1)) { found = false; break; }

        // Set this position as previous for next iteration.
        cp_pos = cn_pos;
    }
    return found;
}


bool ExistWord(const std::string &phrase, const std::string &word)
{
    std::string low_case_phrase = phrase;
    std::transform(std::begin(low_case_phrase), std::end(low_case_phrase),
                   std::begin(low_case_phrase), [](unsigned char c){ return std::tolower(c); });

    std::string low_case_word = word;
    std::transform(std::begin(low_case_word), std::end(low_case_word),
                   std::begin(low_case_word), [](unsigned char c){ return std::tolower(c); });

    return ALGORITHM::ExistExactWord(low_case_phrase, low_case_word);
}


} // End of namespace ALGORITHM
} // End of namespace ELECTRA