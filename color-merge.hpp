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
    std::string output_prefix = "";
} parameters_t;


void parse_arguments(int argc, char **argv, parameters_t & params);

class SDIter {
    sdsl::sd_vector<> *m_s;
    size_t lowndx = 0;
    size_t highndx = 0; // points to next chunk to read
    size_t cur_high = 0;
    size_t nextlow = 0;
    size_t nexthigh = 0;
    size_t next_cur_high = 0;
    size_t num_advanced = 0;
    bool peeked = false;
public:
    SDIter(sdsl::sd_vector<> *a_s) {m_s = a_s;};
    signed long long peek() {
        if (num_advanced >= m_s->low.size()) {
            return -1;
        }
        nextlow = lowndx;
        nexthigh = highndx;
        next_cur_high = cur_high;
        
        size_t high_increment = 0;
        while (m_s->high[nexthigh] != 1) {
            high_increment++;
            nexthigh++;
        }
        nexthigh++; // advance to next one
        next_cur_high += high_increment;
        size_t low = m_s->low[nextlow];
        // for (int i = 1; i <= m_s->wl; i++) {
        //     low |= (m_s->low[nextlow] << (m_s->wl - i));
        // }
        nextlow++;

        peeked = true;
        size_t one_loc = (next_cur_high << m_s->wl) | low; // '1' location
        assert((*m_s)[one_loc]);
        return one_loc;
    };
    void advance() {
            if (!peeked) { peek();}
            lowndx = nextlow;
            highndx = nexthigh;
            cur_high = next_cur_high;
            peeked = false;
            num_advanced++;
    }
};

#endif
