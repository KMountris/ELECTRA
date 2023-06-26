/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#ifndef ELECTRA_APPS_TOOLS_PARSER_TPP_
#define ELECTRA_APPS_TOOLS_PARSER_TPP_

#include "parser.hpp"

namespace APP_ELECTRA
{


template <class T>
T Parser::GetValue(const std::string& attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            std::string error_message = "Required attribute is missing from JSON file: " + attribute;
            throw std::runtime_error(ELECTRA::Logger::Error(error_message));
        }
    }

    T return_value;
    if (current.is_array() && current.size()==1) return_value = current[0];
    else return_value = current;
    return return_value;
}

} // end of namespace APP_ELECTRA


#endif //ELECTRA_APPS_TOOLS_PARSER_TPP_