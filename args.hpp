// args.hpp - tiny command-line helper: --name value / --flag
#pragma once
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

struct Args {
    int argc; char** argv;
    Args(int c, char** v) : argc(c), argv(v) {}
    const char* find(const char* name) const {
        for (int i = 1; i < argc; ++i)
            if (!std::strcmp(argv[i], name)) return i + 1 < argc ? argv[i + 1] : "";
        return nullptr;
    }
    bool has(const char* name) const { return find(name) != nullptr; }
    std::string str(const char* name, const std::string& def = "") const {
        const char* s = find(name); return s ? s : def;
    }
    long num(const char* name, long def) const { const char* s = find(name); return s ? std::atol(s) : def; }
    double real(const char* name, double def) const { const char* s = find(name); return s ? std::atof(s) : def; }
    std::string need(const char* name) const {
        const char* s = find(name);
        if (!s || !*s) throw std::runtime_error(std::string("missing required option ") + name);
        return s;
    }
};
