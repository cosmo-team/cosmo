#ifndef PACK_COLOR_HPP
#define PACK_COLOR_HPP

int getMilliCount();
int getMilliSpan(int nTimeStart);
typedef struct p
{
    std::string plan_filename = "";
    std::string matrix1_filename = "";
    int num_colors1;
    std::string matrix2_filename = "";
    int num_colors2;
    unsigned long long m1;
    unsigned long long n1;
    unsigned long long m2;
    unsigned long long n2;
    std::string output_prefix = "";
} parameters_t;


void parse_arguments(int argc, char **argv, parameters_t & params);


#endif
