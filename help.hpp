#ifndef HPP_HELP
#define HPP_HELP

#include <string>

const std::string g_strHelp = R"(
usage: ./profanity [OPTIONS]

  Basic modes:
    --benchmark             Run without any scoring, a benchmark.
    --zeros                 Score on zeros anywhere in hash.
    --letters               Score on letters anywhere in hash.
    --numbers               Score on numbers anywhere in hash.

  Modes with arguments:
    --leading <single hex>  Score on hashes leading with given hex character.
    --matching <hex string> Score on hashes matching given hex string.

  Advanced modes:
    --leading-range         Scores on hashes leading with characters within
                            given range.
    --range                 Scores on hashes having characters within given
                            range anywhere.

  Range:
    -m, --min <0-15>        Set range minimum (inclusive), 0 is '0' 15 is 'f'.
    -M, --max <0-15>        Set range maximum (inclusive), 0 is '0' 15 is 'f'.

  Misc:
    -s, --skip <index>      Skip device given by index.
    -w, --work <size>       Set OpenCL local work size. [default = 64]
    -W, --workmax <size>    Set OpenCL maximum work size. [default = 32768]

  Examples:
    ./profanity --leading f
    ./profanity --matching dead
    ./profanity --matching badXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXbad
    ./profanity --leading-range -m 0 -M 1
    ./profanity --leading-range -m 10 -M 12
    ./profanity --range -m 0 -M 1

  About:
    profanity is a vanity address generator for Ethereum that utilizes
    computing power from GPUs by using OpenCL.
    
    Author: Johan Gustafsson <profanity@johgu.se>
	Beer donations: 0xbaddec0de82321894ecb1cd4c6d41b831ea15be2)";

#endif /* HPP_HELP */