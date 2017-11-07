#ifndef HPP_CONSTANTS
#define HPP_CONSTANTS

const size_t g_worksizes[] = {
	1,
	255,
	255 * 255,
	255 * 255 * 255
};

#define PROFANITY_PASSES 3
#define PROFANITY_SIZE (g_worksizes[PROFANITY_PASSES])
#define PROFANITY_DEBUG true

#endif /* HPP_CONSTANTS */