/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * ALL RIGHTS RESERVED
 *
 */

#include "parser.hpp"

namespace APP_ELECTRA
{


Parser::Parser(const std::string &json_filename) : json_(), parent_path_("")
{
    std::filesystem::path filepath(json_filename);

    // Check json_filename for *.json extension.
    if (filepath.extension() != ".json") {
        throw std::invalid_argument(ELECTRA::Logger::Error("Could not load configuration file. Expected JSON format (*.json)"));
    }

    // Store parent path of *.json file.
    this->parent_path_ = filepath.parent_path();

    // Load the json file in a stream.
    std::ifstream file_stream(json_filename);

    if (file_stream.is_open()) {
        // Copy the file stream in a string.
        std::string file_str((std::istreambuf_iterator<char>(file_stream)), std::istreambuf_iterator<char>());

        // Parse the file string in the json object.
        this->json_ = nlohmann::json::parse(this->RemoveCommentLines(file_str));

    } else {
        std::string error_msg = "Could not parse ELECTRA application JSON input file. Check file path.";
        throw std::invalid_argument(ELECTRA::Logger::Error(error_msg));
    }

}


Parser::~Parser()
{}


std::string Parser::RemoveCommentLines(const std::string &str)
{
	// Flags indicating single line or multpile line comments have started or not.
	bool sl_cmt = false;
	bool ml_cmt = false;

	// Traverse the input string.
    std::string result;
	for (std::size_t i=0; i != str.length(); ++i) {
		// If single line comment flag is on, then check for its end.
		if (sl_cmt == true && str[i] == '\n') { sl_cmt = false; }

		// If multiple line comment is on, then check for its end.
		else if (ml_cmt == true && str[i] == '*' && str[i+1] == '/') { ml_cmt = false; i++; }

		// If this character is in a comment, ignore it.
		else if (sl_cmt || ml_cmt) { continue; }

		// Check for beginning of comments and set the appropriate flag.
		else if (str[i] == '/' && str[i+1] == '/') { sl_cmt = true; i++; }
		else if (str[i] == '/' && str[i+1] == '*') { ml_cmt = true, i++; }

		// If current character is a non-comment character, append it to result string.
		else result += str[i];
	}

	return result;
}


nlohmann::json Parser::GetObject(const std::string &attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            std::string error_message = attribute + " attribute was not found in the JSON file.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_message));
        }
    }

    return current;
}


bool Parser::HasAttribute(const std::string& attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            return false;
        }
    }

    return true;
}


bool Parser::IsSingleValue(const std::string& attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            std::string error_message = attribute + " attribute was not found in the JSON file.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_message));
        }
    }

    if (current.size() == 1) return true;
    return false;
}


bool Parser::IsArray(const std::string& attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            std::string error_message = attribute + " attribute was not found in the JSON file.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_message));
        }
    }

    if (current.size() > 1 && current[0].size() == 1) return true;
    return false;
}


bool Parser::IsMultiArray(const std::string& attribute) const
{
    std::stringstream test(attribute);
    std::string key_name;
    auto current = this->json_;

    while(std::getline(test, key_name, '.')) {
        nlohmann::json::const_iterator it = current.find(key_name);
        if (it != current.end()) {
            current = it.value();
        } else {
            std::string error_message = attribute + " attribute was not found in the JSON file.";
            throw std::runtime_error(ELECTRA::Logger::Error(error_message));
        }
    }

    if (current.size() > 1 && current[0].size() > 1) return true;
    return false;
}


std::string Parser::ResolvePath(const std::string &filename) const
{
    std::string resolved_filename = filename;
    std::filesystem::path filepath(resolved_filename);
    std::string parent_path_string{filepath.parent_path()};
    if (!filepath.has_root_path() || parent_path_string.size() == 1) {
        if (resolved_filename[0] == '/' || resolved_filename[0] == '\\') {
            resolved_filename = this->parent_path_ + resolved_filename;
        } else {
            resolved_filename = this->parent_path_ + "/" + resolved_filename;
        }
    }

    return resolved_filename;
}


} // end of namespace APP_ELECTRA