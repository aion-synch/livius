#ifndef CORENC_COLORS_H
#define CORENC_COLORS_H
#include <string>
#include <iostream>
namespace corenc
{
    namespace color
    {
        const std::string ESCAPE = "\u001b[0m";
        // 8-bit colors
        const std::string BLACK = "\u001b[30m";
        const std::string RED = "\u001b[31m";
        const std::string GREEN = "\u001b[32m";
        const std::string YELLOW = "\u001b[33m";
        const std::string BLUE = "\u001b[34m";
        const std::string MAGENTA = "\u001b[35m";
        const std::string CYAN = "\u001b[36m";
        const std::string WHITE = "\u001b[37m";
        const std::string PURPLE = "\e[1;35m";
        // 16-bit colors
        const std::string BBLACK = "\u001b[30;1m";
        const std::string BRED = "\u001b[31;1m";
        const std::string BGREEN = "\u001b[32;1m";
        const std::string BYELLOW = "\u001b[33;1m";
        const std::string BBLUE = "\u001b[34;1m";
        const std::string BMAGENTA = "\u001b[35;1m";
        const std::string BCYAN = "\u001b[36;1m";
        const std::string BWHITE = "\u001b[37;1m";
        static void color_output(const std::string& text, const std::string& col)
        {
            std::cout << col << text << ESCAPE << std::endl;
        }
        static void color_output(std::ostream& os, const std::string& text, const std::string& col)
        {
            os << col << text << ESCAPE << std::endl;
        }
    }
}
#endif // CORENC_COLORS_H
