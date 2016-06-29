#ifndef PACK_COLOR_HPP
#define PACK_COLOR_HPP

int getMilliCount();
int getMilliSpan(int nTimeStart);
typedef struct p
{
    std::string input_filename = "";
    int num_colors;
    unsigned long long m;
    unsigned long long n;
    std::string output_prefix = "";
} parameters_t;


void parse_arguments(int argc, char **argv, parameters_t & params);


#endif
